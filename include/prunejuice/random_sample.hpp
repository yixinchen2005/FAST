#pragma once

#include <boost/bimap.hpp>
#include <havoqgt/mpi.hpp>
#include <vector>
#include <set>
#include <algorithm>
#include "query_cache.hpp"
// timer
#include <signal.h>    /* union sigval / struct sigevent */
#include <stdio.h>    /* printf */
#include <string.h>    /* memset */
#include <unistd.h> /* sleep */
#include <time.h>

#include <ctime>
#include <cstdlib>

#define NULL_NODE 0xffffffffffffffff

#define MAX_WALK_SIZE 50

// #define ENABLE_CACHE

using namespace havoqgt;


namespace prunejuice
{
    typedef uint64_t Count;
    typedef uint64_t Vertex;
    typedef uint64_t Edge;
    typedef uint64_t VertexData;
    typedef uint64_t VertexClass;
    typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;
    typedef graph_type::vertex_locator vertex_locator;

    template <typename Visitor>
    class Sample_queue_tds{
    public:
        Sample_queue_tds() {}

        bool push(Visitor const &element)
        {
            // std::cout<<"push"<<std::endl;
            data.push_back(element);
            return true;
        }

        void pop()
        {
            //data.pop_back();
            data.pop_front();
        }

        Visitor const &top()
        {
            //return data.back();
            return data.front();
        }

        size_t size() const
        {
            return data.size();
            ;
        }

        bool empty() const
        {
            return data.empty();
        }

        void clear()
        {
            data.clear();
        }

    protected:
        //std::vector<Visitor> data;
        std::deque<Visitor> data;
    };

    /**
 * @brief visitor
 */
    template <typename Graph>
    class Sample_visitor
    {
    public:
        template <typename Array>
        Sample_visitor(vertex_locator vertex_, Vertex v, Array& walk_, int& walk_size_, int msg_, bool whether_add_path){
            vertex = vertex_;
            msg = msg_;
            walk_size = walk_size_;
            std::copy(std::begin(walk_), std::end(walk_),std::begin(walked_path));
            if(whether_add_path){
                walked_path[walk_size_]=v;
                walk_size = walk_size_+1;
            }
        }

        // template <typename Array>
        // Sample_visitor(vertex_locator vertex_, Vertex v, Array& walk_, int& walk_size_, Array& label_path_){
        //     vertex = vertex_;
        //     std::copy(std::begin(walk_), std::end(walk_),std::begin(walked_path));
        //     std::copy(std::begin(label_path_), std::end(label_path_),std::begin(label_path));
        //     walked_path[walk_size_]=v;
        //     walk_size = walk_size_+1;
        // }

        Sample_visitor(vertex_locator vertex_){
            vertex = vertex_;
            walk_size = 0;
            msg = 0;
        }
        Sample_visitor(){
            walk_size = 0;
            msg = 0;
        }

        template <typename AlgData>
        bool pre_visit(AlgData &alg_data) const{
            return true;
        }

        template <typename VisitorQueueHandle, typename AlgData>
        bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data){
            walked_path[walk_size] = g.locator_to_label(vertex);
            ++ walk_size;
            return visit(g, vis_queue, alg_data);
        }

        template <typename VisitorQueueHandle, typename AlgData>
        bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
        {
            double time_start,time_end;
            if(msg == 0){
                // randomly select a neighbor
                if(walk_size >= MAX_WALK_SIZE){
                    for(int i=0;i<walk_size;++i){
                        Vertex v = walked_path[i];
                        vertex_locator target_locator = g.label_to_locator(v);
                        Sample_visitor new_visitor(target_locator, v, walked_path, walk_size, 1, false);
                        vis_queue->queue_visitor(new_visitor);
                    }
                    return true;
                }
                std::unordered_set<Vertex> can;
                std::unordered_set<Vertex> path;
                for(int i=0;i<walk_size;++i){
                 path.insert(walked_path[i]);
                }
                for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
                    vertex_locator neighbor = eitr.target();
                    if(path.find(g.locator_to_label(neighbor)) != path.end()){
                     continue;
                    }
                    can.insert(g.locator_to_label(neighbor));
                }
                
                if(can.size() == 0){
                    // return to a path randomly
                    int select = rand()%path.size();
                    if(path.size()){
                        std::cout<<"ERROR: path_size is 0"<<std::endl;
                        exit(0);
                    }
                    for(auto v : path){
                        if(select == 0){
                            // send to the origin
                            vertex_locator target_locator = g.label_to_locator(v);
                            Sample_visitor new_visitor(target_locator, v, walked_path, walk_size, 0, false);
                            vis_queue->queue_visitor(new_visitor);
                            return true;
                        }
                        -- select;
                    }
                }
                int select = rand()%can.size();
                for(auto v: can){
                    if(select == 0){
                        vertex_locator target_locator = g.label_to_locator(v);
                        Sample_visitor new_visitor(target_locator, v, walked_path, walk_size, 0, true);
                        vis_queue->queue_visitor(new_visitor);
                        return true;
                    }
                    -- select;
                }
            }else if(msg == 1){
                // output the graph
                auto & vertex_label = std::get<1>(alg_data);
                auto & label_map = std::get<2>(alg_data);
                auto & edge_list = std::get<3>(alg_data);
// 
                std::unordered_set<Vertex> path;
                for(int i=0;i<walk_size;++i){
                    path.insert(walked_path[i]);
                }
                label_map.insert({g.locator_to_label(vertex), vertex_label[vertex]});
                for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
                    vertex_locator neighbor = eitr.target();
                    if(path.find(g.locator_to_label(neighbor)) == path.end()){
                     continue;
                    }
                    edge_list.push_back({g.locator_to_label(vertex), g.locator_to_label(neighbor)});
                }
            }
            
            return true;
        }
        
       //  template<typename AlgData>
       //  bool check_label(AlgData& alg_data){
       //      if(walk_size>1){
       //          if(label_path[walk_size-1] == label_path[walk_size-2]){
       //              return false;
       //          }
       //          if(label_path[walk_size-1]!=2 && label_path[walk_size-1]!=3 && label_path[walk_size-1]!=4 && label_path[walk_size-1]!=5){
       //              return false;
       //          }
       //          std::unordered_set<Vertex> labels;
       //          for(int i=0;i<walk_size-1;++i){
       //              labels.insert(label_path[i]);
       //          }
       //          if(labels.size()<3 && walk_size-1>3){
       //              return false;
       //          }
       //      }else{
       //          if(label_path[walk_size-1]!=2 && label_path[walk_size-1]!=3 && label_path[walk_size-1]!=4 && label_path[walk_size-1]!=5){
       //              return false;
       //          }
       //      }
       //      return true;
       //  }

       //  template <typename VisitorQueueHandle, typename AlgData>
       //  bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
       //  {
       //      // preprocess the label
       //      auto label = std::get<1>(alg_data)[vertex];
       //      label_path[walk_size-1] = label;
       //      if(!check_label(alg_data)){
       //          return false;
       //      }

       //      if(walk_size == MAX_WALK_SIZE){
       //          // check the connection with the endpoint
       //          for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
       //              vertex_locator neighbor = eitr.target();
       //              if(g.locator_to_label(neighbor) == walked_path[0]){
       //                  for(int i=0;i<walk_size;++i){
       //                      std::get<4>(alg_data)<<walked_path[i]<<":"<<label_path[i]<<" ";
       //                  }
       //                  std::get<4>(alg_data)<<std::endl;
       //                  std::cout<<"find!"<<std::endl;
       //                  MPI_Abort(MPI_COMM_WORLD, -1);
       //                  return false;
       //              }
       //          }
       //      }else{
       //          for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
       //              vertex_locator neighbor = eitr.target();
       //              Sample_visitor new_visitor(neighbor, g.locator_to_label(neighbor), walked_path, walk_size, label_path);
       //              vis_queue->queue_visitor(new_visitor);
       //          }
       //      }
       //      return true;
       //  }


        friend inline bool operator>(const Sample_visitor &v1, const Sample_visitor &v2)
        {
            if (v1.walk_size > v2.walk_size){
                return false;
            }
            return true;
        }

        vertex_locator vertex;
        std::array<Vertex, MAX_WALK_SIZE> walked_path; // 16 is the maximal query size
        std::array<Vertex, MAX_WALK_SIZE> label_path; // 16 is the maximal query size
        int msg; // 0-has neighbor in the same partition otherwise 1
        int walk_size;
    };

    /**
 * @brief start enumeration 
 */
    template <typename DataGraph, typename VertexMetadata, typename LabelMap, typename EdgeList, typename File>
    void sample_graph(DataGraph *graph, VertexMetadata vertex_label, int local_node_count, LabelMap& label_map, EdgeList& edge_list, File& result)
    {
        auto end = false;
        auto alg_data = std::forward_as_tuple(graph, vertex_label, label_map, edge_list, result, end);
        int mpi_rank = havoqgt::comm_world().rank();
        


        typedef Sample_visitor<DataGraph> visitor_type;
        auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(graph, alg_data);
        // auto vq = create_visitor_queue<visitor_type, Sample_queue_tds>(graph, alg_data);
        MPI_Barrier(MPI_COMM_WORLD);

        std::vector<vertex_locator> start_visitors;
        // find the label of the first query vertex
        if(mpi_rank == 0){
            std::srand((unsigned)std::time(NULL));
            int select = rand()%local_node_count;
            int idx=0;
            for (auto vitr = graph->vertices_begin(); vitr != graph->vertices_end();++vitr, ++idx){
                if(idx == select){
                    start_visitors.push_back(*vitr);
                    break;
                }
            }
        }
        vq.init_visitor_traversal(start_visitors);
    }
} // namespace prunejuice
