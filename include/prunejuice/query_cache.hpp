#pragma once

#include <boost/bimap.hpp>
#include <havoqgt/mpi.hpp>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_set>
// timer
#include <signal.h>    /* union sigval / struct sigevent */
#include <stdio.h>    /* printf */
#include <string.h>    /* memset */
#include <unistd.h> /* sleep */
#include <time.h>

#define CACHE_SIZE 20000 // cache size
#define UPDATE_BATCH_SIZE 100 // cache line size, unit: number of vertices
#define UNHIT_LIMIT 10000 // after a missing cache line unhitted over 10000 times, the cache is force to updated
#define PREFECT_RATIO 2 // 1/2 of the whole content is loaded in advance

typedef uint64_t Vertex;
typedef uint64_t Edge;

class Query_cache
{
private:
    // std::vector<Vertex> offset;
    // std::vector<Vertex> neighbor_list;
    Vertex* offset;
    Vertex* neighbor_list;
    MPI_Win win;
    MPI_Win win_neighor;
    int mpi_size;
    int mpi_rank;
    std::unordered_map<Vertex, std::vector<Vertex>> cache;
    Vertex vertex_count;
public:
    // metrics
    Vertex hit_count, total_count, unhit_count;
    std::unordered_set<Vertex> request;

    template <typename VertexActive, typename EdgeActive, typename Graph>
    Query_cache(VertexActive& active_vertex, EdgeActive& active_edge, Graph& graph, int mpi_rank, int mpi_size, long long unsigned vertex_count){
        Vertex offset_index = 0;
        Vertex edge_count = 0;
        hit_count = 0;
        total_count = 0;
        unhit_count = 0;
        this->mpi_size = mpi_size;
        this->mpi_rank = mpi_rank;
        this->vertex_count = vertex_count;

        for(int i=mpi_rank;i<vertex_count;i+=mpi_size){
            auto v_locator = graph->label_to_locator(i);
            edge_count += active_edge[v_locator].size();
        }

        MPI_Win_allocate((MPI_Aint)((vertex_count/mpi_size+1) * sizeof(Vertex)), sizeof(Vertex),
                 MPI_INFO_NULL, MPI_COMM_WORLD, &offset, &win);
        MPI_Win_allocate((MPI_Aint)((edge_count+1) * sizeof(Vertex)), sizeof(Vertex),
                 MPI_INFO_NULL, MPI_COMM_WORLD, &neighbor_list, &win_neighor);
        Vertex acc = 0;
        int idx = 0;
        for(int i=mpi_rank;i<vertex_count;i+=mpi_size, ++idx){
            offset[idx] = offset_index;
            auto v_locator = graph->label_to_locator(i);
            for (auto &n : active_edge[v_locator])
            {
                neighbor_list[acc] = n.first;
                ++ acc;
            }
            // reorder
            for(int i=offset_index; i<offset_index+active_edge[v_locator].size(); ++i){
                for(int j=i+1;j<offset_index+active_edge[v_locator].size(); ++j){
                    if(neighbor_list[i]>neighbor_list[j]){
                        std::swap(neighbor_list[i], neighbor_list[j]);
                    }
                }
            }
            offset_index += active_edge[v_locator].size();
        }
        offset[idx] = offset_index;
        
        MPI_Barrier(MPI_COMM_WORLD);
        pre_fetch();
    }

    void get_local_neighbor(Vertex v){
        // std::cout<<"-->"<<v/mpi_size<<"-:-"<<offset[v/mpi_size]<<std::endl;
        for(int i=offset[v/mpi_size];i<offset[v/mpi_size+1]; ++i){
            std::cout<<mpi_rank<<":"<<neighbor_list[i]<<std::endl;
        }
    }

    std::pair<Vertex*, Vertex> get_neighbor(Vertex& v){
        if(v%mpi_size == mpi_rank){
            Vertex start_offset = v/mpi_size;
            return {neighbor_list+offset[start_offset], offset[start_offset+1]-offset[start_offset]};
        }else{
            auto itf = cache.find(v);
            ++ total_count;
            if(itf != cache.end()){
                ++ hit_count;
                // std::cout<<"hit"<<std::endl;
                return {&(itf->second[0]), itf->second.size()};   
            }else{
                request.insert(v);
                if((request.size())%UPDATE_BATCH_SIZE == 0 || unhit_count > UNHIT_LIMIT){
                    // std::cout<<"updating"<<std::endl;
                    get_remote_neighbor_list(request);
                    request.clear();
                    itf = cache.find(v);
                    unhit_count = 0;
                    return {&(itf->second[0]), itf->second.size()}; 
                }
                ++ unhit_count;
                return {NULL, -1};
            }
            
        }
    }

    ~Query_cache(){
        // MPI_Win_wait(win);
        // MPI_Win_wait(win_neighor);
        MPI_Win_free(&win);
        MPI_Win_free(&win_neighor);
    }

    void pre_fetch(){
        Vertex idx = 0;
        std::unordered_set<Vertex> prefetch;
        for(int i=mpi_rank; i<vertex_count; i+=mpi_size, ++idx){
            Vertex n_size = offset[idx+1]-offset[idx];
            if(n_size > 0){
                for(int j=offset[idx]; j<offset[idx+1]; j++){
                    prefetch.insert(neighbor_list[j]);
                    if(prefetch.size() >= CACHE_SIZE/PREFECT_RATIO){
                        goto out;
                    }
                }
            }
        }
        out:
        get_remote_neighbor_list(prefetch);
    }

    void get_remote_neighor(Vertex v){
        int target_rank = v%mpi_size;
        // MPI_Win_start(MPI_COMM_WORLD, 0, win);
        MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0, win);
        Vertex buf[2];
        MPI_Get(&buf, 2, MPI_UNSIGNED_LONG_LONG, 0, 0, 2, MPI_UNSIGNED_LONG_LONG, win);
        // MPI_Win_complete(win);
        MPI_Win_unlock(target_rank, win);
        // MPI_Win_start(MPI_COMM_WORLD, 0, win_neighor);
        
        Vertex neighbor_size = buf[1]-buf[0];
        std::vector<Vertex> b;
        if(neighbor_size > 0){
            b.resize(neighbor_size);
            MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0, win_neighor);
            MPI_Get(&b[0], neighbor_size, MPI_INT64_T, target_rank, buf[0], neighbor_size, MPI_INT64_T, win_neighor);
            MPI_Win_unlock(target_rank, win_neighor);
            for(int i=0;i<neighbor_size;i++){
                std::cout<<b[i]<<std::endl;
            }
        }
    }

    void get_remote_neighbor_list(std::unordered_set<Vertex>& source_list){
        std::vector<Vertex> buf_size;
        buf_size.resize(source_list.size()*2);
        Vertex idx = 0;
        Vertex total_size = 0;
        for(auto v:source_list){
            int target_rank = v%mpi_size;
            MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0, win);
            // std::cout<<v<<"b:"<<buf_size[idx*2]<<":"<<buf_size[idx*2+1]<<std::endl;
            MPI_Get(&buf_size[idx*2], 2, MPI_INT64_T, target_rank, (v/mpi_size), 2, MPI_INT64_T, win);
            MPI_Win_unlock(target_rank, win);
            total_size += buf_size[idx*2+1]-buf_size[idx*2];
            // std::cout<<v<<":"<<buf_size[idx*2]<<":"<<buf_size[idx*2+1]<<std::endl;
            ++idx;
        }
        
        std::vector<Vertex> recv_buf;
        recv_buf.resize(total_size);
        idx = 0;
        Vertex acc_size = 0;
        for(auto v:source_list){
            Vertex neighbor_size = buf_size[idx*2+1]-buf_size[idx*2];
            if(neighbor_size>0){
                int target_rank = v%mpi_size;
                MPI_Win_lock(MPI_LOCK_SHARED, target_rank, 0, win_neighor);
                MPI_Get(&recv_buf[acc_size], neighbor_size, MPI_INT64_T, target_rank, buf_size[idx*2], neighbor_size, MPI_INT64_T, win_neighor);
                MPI_Win_unlock(target_rank, win_neighor);
                buf_size[idx*2] = acc_size;
                acc_size += neighbor_size;
                buf_size[idx*2+1] = acc_size;
            }else{
                buf_size[idx*2] = acc_size;
                buf_size[idx*2+1] = acc_size;
            }
            ++ idx;
        }

        add_cache_list(source_list, buf_size, recv_buf);
    }

    void add_cache_list(Vertex& v, std::vector<Vertex>& neighbor){
        if(cache.size()+1>CACHE_SIZE){
            replace_policy(1);
        }
        cache.insert({v, neighbor});
    }

    void add_cache_list(std::unordered_set<Vertex>& source_list, std::vector<Vertex>& offset_list, std::vector<Vertex>& neighbors){
        if(cache.size()+source_list.size()>CACHE_SIZE){
            replace_policy(cache.size()+source_list.size()-CACHE_SIZE);
        }
        Vertex idx = 0;
        for(auto v:source_list){
            std::vector<Vertex> t(std::next(neighbors.begin(), offset_list[idx*2]), std::next(neighbors.begin(), offset_list[idx*2+1]));
            cache.insert({v, t});
            ++idx;
        }
    }

    void replace_policy(Vertex evict_size){
        // random replace
        if(evict_size > cache.size())
            evict_size = cache.size();
        cache.erase(cache.begin(), std::next(cache.begin(), evict_size));
    }

    void print_cache(std::vector<Vertex>& source_list){
        for(auto v:source_list){
            auto itf = cache.find(v);
            std::cout<<v<<":{";
            if(itf!=cache.end()){
                for(auto n:itf->second){
                    std::cout<<n<<",";
                }
            }
            std::cout<<"}"<<std::endl;
        }
    }

    void test_query(Vertex v){
        auto p = get_neighbor(v);
        if(p.second == -1){
            std::cout<<v<<":unhit"<<std::endl;
        }else{
            std::cout<<v<<":{";
            for(int i=0;i<p.second;++i){
                std::cout<<p.first[i]<<",";
            }
            std::cout<<"}"<<std::endl;
        }
    }

};

