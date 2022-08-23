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

#include "config.hpp"

#define NULL_NODE 0xffffffffffffffff

#define MAX_QUERY_SIZE 8

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
    static constexpr size_t max_bit_vector_size = 16; 
    typedef std::bitset<max_bit_vector_size> BitSet;

    template <typename Visitor>
    class Enumeration_queue_tds{
    public:
        Enumeration_queue_tds() {}

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

    timer_t id;
    void set_stop_flag(union sigval val)
    {
        // printf("%d\n",val.sival_int);
        // printf("%p\n",val.sival_ptr);
        *((bool*)val.sival_ptr) = true;
        timer_delete(id);
    }

    void register_timer(bool* stop_flag){
        struct timespec spec;
        struct sigevent ent;
        struct itimerspec value;
        struct itimerspec get_val;

        /* Init */
        memset(&ent, 0x00, sizeof(struct sigevent));
        memset(&get_val, 0x00, sizeof(struct itimerspec));

        int test_val = 0;
        /* create a timer */
        ent.sigev_notify = SIGEV_THREAD;
        ent.sigev_notify_function = set_stop_flag;
        ent.sigev_value.sival_ptr = stop_flag;
        // printf("create timer\n");
        timer_create(CLOCK_MONOTONIC, &ent, &id);

        /* start a timer */
        value.it_value.tv_sec = ENUMERATION_TIME;
        value.it_value.tv_nsec = 0;
        value.it_interval.tv_sec = 0;
        value.it_interval.tv_nsec = 0;
        // printf("start timer\n");
        timer_settime(id, 0, &value, NULL);
    }


    /**
 * @brief visitor
 */
    template <typename Graph>
    class Enumeration_visitor
    {
    public:
#ifdef ENABLE_CACHE
        template <typename Array, typename ArrayConn>
        Enumeration_visitor(vertex_locator vertex_, Array& partial_match_, Vertex v, int partial_size_, int msg_, ArrayConn& Conn_check_)
        :vertex(vertex_), msg(msg_){
            parent_mpi_rank = havoqgt::comm_world().rank();
            std::copy(std::begin(partial_match_), std::end(partial_match_),std::begin(partial_match));
            std::copy(std::begin(Conn_check_), std::end(Conn_check_),std::begin(Conn_check));
            if(msg_ == 0){
                partial_match[partial_size_]=v;
                partial_size = partial_size_+1;
            }else{
                partial_size = partial_size_;
            }
        }

        Enumeration_visitor(vertex_locator vertex_){
            vertex = vertex_;
            partial_size = 0;
            msg = 0;
            memset(&Conn_check[0], MAX_QUERY_SIZE, sizeof(bool)*MAX_QUERY_SIZE);
        }

        Enumeration_visitor(){
            partial_size = 0;
            msg = 0;
            memset(&Conn_check[0], MAX_QUERY_SIZE, sizeof(bool)*MAX_QUERY_SIZE);
        }

#else
        template <typename Array>
        Enumeration_visitor(vertex_locator vertex_, Array& partial_match_, Vertex v, int partial_size_, int msg_)
        :vertex(vertex_), msg(msg_){
            parent_mpi_rank = havoqgt::comm_world().rank();
            std::copy(std::begin(partial_match_), std::end(partial_match_),std::begin(partial_match));
            if(msg_ == 0){
                partial_match[partial_size_]=v;
                partial_size = partial_size_+1;
            }else{
                partial_size = partial_size_;
            }
        }

        Enumeration_visitor(vertex_locator vertex_){
            vertex = vertex_;
            partial_size = 0;
            msg = 0;
        }
        Enumeration_visitor(){
            partial_size = 0;
            msg = 0;
        }
#endif

        template <typename AlgData>
        bool pre_visit(AlgData &alg_data) const{
            if(std::get<14>(alg_data)){
                return false;
            }
            return true;
        }

        template <typename VisitorQueueHandle, typename AlgData>
        bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data){
            partial_match[partial_size] = g.locator_to_label(vertex);
            ++ partial_size;
            parent_mpi_rank = havoqgt::comm_world().rank();
            // for(int i=0;i<partial_size;++i){
            //     std::cout<<partial_match[i]<<" ";
            // }
            // std::cout<<std::endl;
            return visit(g, vis_queue, alg_data);
        }

        inline bool check(Vertex v){
            for(int i=0;i<partial_size;++i){
                if(partial_match[i] == v){
                    return true;
                }
            }
            return false;
        }

        inline bool check_before(Vertex v){
            for(int i=0;i<partial_size-1;++i){
                if(partial_match[i] == v){
                    return true;
                }
            }
            return false;
        }

        template <typename VisitorQueueHandle, typename AlgData>
        bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
        {
            double time_start,time_end;
            // time_start = MPI_Wtime();
            auto& pattern_graph = std::get<0>(alg_data);
            auto& ordering = std::get<1>(alg_data);
            auto& vertex_metadata = std::get<2>(alg_data);
            auto& vertex_state_map = std::get<3>(alg_data);
            auto& template_vertices = std::get<4>(alg_data);
            auto& plan = std::get<5>(alg_data);
            auto& active_vertex_map = std::get<6>(alg_data);
            auto& active_edge_map = std::get<7>(alg_data);
            auto& active_edge = std::get<12>(alg_data);
            auto& time_col = std::get<13>(alg_data);
            auto & stop = std::get<14>(alg_data);
#ifdef ENABLE_CACHE
            auto & query_cache = std::get<16>(alg_data);
#endif
            // generates the candidates  
            if(stop){
                vis_queue->clear();
                return false;
            } 
            Vertex extending_vertex = ordering[partial_size];
            vector<Vertex> candidates;
            int index = 0;
            int mpi_rank = havoqgt::comm_world().rank();
            int next_msg = 1;
            
            if(msg == 0){
                std::get<15>(alg_data)[0] ++;
                // active vertex pruning
                if (active_vertex_map.find(g.locator_to_label(vertex)) == active_vertex_map.end()) {
                    return false;
                }

                // candidate pruning
                auto offset_map = std::get<9>(alg_data);
                BitSet query_candidates(template_vertices[vertex]);
                Vertex current_query_vertex = ordering[partial_size-1];
                if(!query_candidates[current_query_vertex]){
                    return false;
                }

                // checking connection
#ifdef ENABLE_CACHE
                for(int i=0;i<partial_size-1;++i){
                    if(Conn_check[i] == true){
                        Vertex match_vertex = partial_match[i];
                        if(active_edge_map[vertex].find(match_vertex) == active_edge_map[vertex].end()){
                            std::get<15>(alg_data)[2] ++;
                            return false;
                        }
                    }
                }
#else
                for(auto q_n : plan[partial_size-1]){
                    Vertex match_vertex = partial_match[offset_map[q_n]];
                    if(active_edge_map[vertex].find(match_vertex) == active_edge_map[vertex].end()){
                        std::get<15>(alg_data)[2] ++;
                        return false;
                    }
                }
#endif

                // printing the results
                if(partial_size == ordering.size()){
#ifdef PRINT_RESULT
                    for(int i=0;i<partial_size;++i){
                        std::get<8>(alg_data)<<partial_match[offset_map[i]]<<" ";
                        // std::get<8>(alg_data)<<partial_match[i]<<" ";
                    }
                    std::get<8>(alg_data)<<std::endl;
#endif
                    std::get<11>(alg_data) += 1;
                    return true;
                }

#ifdef ENABLE_CACHE
                memset(&Conn_check[0], 0, sizeof(bool)*partial_size);
                for(auto q_n : plan[partial_size]){
                    Vertex match_vertex = partial_match[offset_map[q_n]];
                    vertex_locator match_locator = g.label_to_locator(match_vertex);
                    auto p = query_cache.get_neighbor(match_vertex);
                    if(p.second == -1){
                        Conn_check[offset_map[q_n]] = true;
                        continue;
                    }
                    if(p.first == NULL){
                        return false;
                    }
                    next_msg = 0;
                    if(index == 0){
                        for(int i=0;i<p.second;++i){
                            candidates.push_back(p.first[i]);
                        }
                    }else{
                        std::vector<Vertex> tmp_candidates;
                        std::vector<Vertex> c1(p.first, p.first+p.second);
                        std::set_intersection(candidates.begin(), candidates.end(), c1.begin(), c1.end(), std::insert_iterator<std::vector<Vertex> >(tmp_candidates, tmp_candidates.begin()));
                        std::swap(tmp_candidates, candidates);
                    }
                    ++index;
                }
#else
                for(auto q_n : plan[partial_size]){
                    Vertex match_vertex = partial_match[offset_map[q_n]];
                    vertex_locator match_locator = g.label_to_locator(match_vertex);
                    if(match_locator.owner() != mpi_rank){
                        continue;
                    }
                    next_msg = 0;
                    if(index == 0){
                        for(auto &neighbor : active_edge_map[match_locator]){
                            candidates.push_back(neighbor.first);
                        }
                    }else{
                        vector<Vertex> tmp_candidates;
                        for(auto v: candidates){
                            auto itf = active_edge_map[match_locator].find(v);
                            if(itf != active_edge_map[match_locator].end()){
                                tmp_candidates.push_back(v);
                            }
                        }
                        swap(candidates, tmp_candidates);
                    }
                    ++index;
                }
#endif

// #ifdef ENABLE_CACHE
//                 memset(&Conn_check[0], 0, sizeof(bool)*partial_size);
// #endif
//                 // time_start = MPI_Wtime();
//                 for(auto q_n : plan[partial_size]){
//                     Vertex match_vertex = partial_match[offset_map[q_n]];
//                     vertex_locator match_locator = g.label_to_locator(match_vertex);
//                     if(match_locator.owner() != mpi_rank){
// #ifdef ENABLE_CACHE
//                         Conn_check[offset_map[q_n]] = true;
// #endif
//                         continue;
//                     }
//                     next_msg = 0;
//                     if(index == 0){
//                         for(auto &neighbor : active_edge_map[match_locator]){
//                             candidates.push_back(neighbor.first);
//                         }
//                     }else{
//                         vector<Vertex> tmp_candidates;
//                         for(auto v: candidates){
//                             auto itf = active_edge_map[match_locator].find(v);
//                             if(itf != active_edge_map[match_locator].end()){
//                                 tmp_candidates.push_back(v);
//                             }
//                         }
//                         swap(candidates, tmp_candidates);
//                     }
//                     // std::get<8>(alg_data)<<"-------------:"<<match_vertex<<std::endl;
//                     // for(auto v:candidates){
//                     //     std::get<8>(alg_data)<<v<<" ";
//                     // }
//                     // std::get<8>(alg_data)<<std::endl;
//                     ++index;
//                 }


                // time_start = MPI_Wtime();
                if(next_msg == 0){
                    // time_start = MPI_Wtime();
                    for(auto c : candidates){
                        if(check(c)){
                            continue;
                        }
                        // if(c<partial_match[partial_size-1]){
                        //     continue;
                        // }
                        vertex_locator c_locator = g.label_to_locator(c);
#ifdef ENABLE_CACHE
                        Enumeration_visitor new_visitor(c_locator, partial_match, c, partial_size, 0, Conn_check);
#else
                        Enumeration_visitor new_visitor(c_locator, partial_match, c, partial_size, 0);
#endif
                        // Enumeration_visitor new_visitor(c_locator, partial_match, c, partial_size, 0);
                        vis_queue->queue_visitor(new_visitor);
                    }
                    // time_end = MPI_Wtime();
                    // time_col[2] += time_end-time_start;
                }else{
                    Vertex q_n = plan[partial_size][0];
                    Vertex match_vertex = partial_match[offset_map[q_n]];
                    vertex_locator match_locator = g.label_to_locator(match_vertex);
#ifdef ENABLE_CACHE
                    Enumeration_visitor new_visitor(match_locator, partial_match, 0, partial_size, 1, Conn_check);
#else
                    Enumeration_visitor new_visitor(match_locator, partial_match, 0, partial_size, 1);
#endif
                    // Enumeration_visitor new_visitor(match_locator, partial_match, 0, partial_size, 1);
                    vis_queue->queue_visitor(new_visitor);
                }
            }else{
                // time_start = MPI_Wtime();
                std::get<15>(alg_data)[1] ++;
                for(auto &neighbor : active_edge_map[vertex]){
                    if(check(neighbor.first)){
                        continue;
                    }
                    vertex_locator match_locator = g.label_to_locator(neighbor.first);
#ifdef ENABLE_CACHE
                    Enumeration_visitor new_visitor(match_locator, partial_match, neighbor.first, partial_size, 0, Conn_check);
#else
                    Enumeration_visitor new_visitor(match_locator, partial_match, neighbor.first, partial_size, 0);
#endif
                    // Enumeration_visitor new_visitor(match_locator, partial_match, neighbor.first, partial_size, 0);
                    vis_queue->queue_visitor(new_visitor);
                }
                // time_end = MPI_Wtime();
                // time_col[3] += time_end-time_start;
            }
            return true;
        }

        friend inline bool operator>(const Enumeration_visitor &v1, const Enumeration_visitor &v2)
        {
            if (v1.partial_size > v2.partial_size){
                return false;
            }
            return true;
        }

        vertex_locator vertex;
        std::array<Vertex, MAX_QUERY_SIZE> partial_match; // 16 is the maximal query size
#ifdef ENABLE_CACHE
        std::array<bool, MAX_QUERY_SIZE> Conn_check;
#endif 
        int msg; // 0-has neighbor in the same partition otherwise 1
        int partial_size;
        int parent_mpi_rank;
    };

    /**
 * @brief start enumeration 
 */
#ifdef ENABLE_CACHE
    template <typename DataGraph, typename PatternGraph, typename VertexMetadata, typename VertexStateMap, typename TemplateVertex, typename VertexActive, typename EdgeActive, typename ActiveMap, typename timeCol, typename LocalCompute, typename QCACHE>
    void pattern_enumeration(DataGraph *graph, PatternGraph* pattern_graph, std::vector<Vertex>& ordering, VertexMetadata& vertex_metadata,
                             VertexStateMap& vertex_state_map, TemplateVertex& template_vertices, VertexActive vertex_active_state_map, EdgeActive vertex_active_edges_map, std::ofstream& result, unsigned long long& result_count, ActiveMap& active_edge_map, timeCol& time_col, LocalCompute& localCompute, QCACHE& query_cache, std::ofstream& pruned_graph)
    
#else
    template <typename DataGraph, typename PatternGraph, typename VertexMetadata, typename VertexStateMap, typename TemplateVertex, typename VertexActive, typename EdgeActive, typename ActiveMap, typename timeCol, typename LocalCompute>
    void pattern_enumeration(DataGraph *graph, PatternGraph* pattern_graph, std::vector<Vertex>& ordering, VertexMetadata& vertex_metadata,
                             VertexStateMap& vertex_state_map, TemplateVertex& template_vertices, VertexActive vertex_active_state_map, EdgeActive vertex_active_edges_map, std::ofstream& result, unsigned long long& result_count, ActiveMap& active_edge_map, timeCol& time_col, LocalCompute& localCompute, std::ofstream& pruned_graph)
    
#endif
    {
        // typedef Enumeration_visitor<> visitor_type;
        // generate query plan
        std::vector<std::vector<Vertex>> plan;
        // key:vertex_id, value:offset in ordering
        std::unordered_map<Vertex, int> offset_map;
        for(int i=0;i<ordering.size();++i){
            offset_map.insert({ordering[i], i});
        }
        bool stop=false;
        
        // for(int i=0;i<ordering.size();++i){
        //     std::cout<<ordering[i]<<":";
        //     for(auto v : pattern_graph->neighbor[ordering[i]]){
        //         std::cout<<v<<" ";
        //     }
        //     std::cout<<std::endl;
        // }
        plan.push_back({});
        for(int i=1;i<ordering.size();++i){
            Vertex current_vertex = ordering[i];
            std::vector<Vertex> tmp;
            plan.push_back(tmp);
            for(int j=0;j<i;++j){
                for(auto n : pattern_graph->neighbor[ordering[i]]){
                    if(n == ordering[j]){
                        plan[i].push_back(n);
                    }
                }
            }
            // std::cout<<"plan-"<<ordering[i]<<":";
            // for(auto v : plan[i]){
            //     std::cout<<v<<" ";
            // }
            // std::cout<<std::endl;
        }
        // printing the pruned graph
        // unsigned long long u_count=0,e_count=0;
        // for(auto& v:vertex_state_map){
        //      pruned_graph<<"v "<<v.first<<" "<<vertex_metadata[graph->label_to_locator(v.first)]<<std::endl;
        // }
        // for(auto& v:vertex_state_map){
        //      u_count ++;
        //      for(auto& neighbor : vertex_active_edges_map[graph->label_to_locator(v.first)]){
        //          pruned_graph<<"e "<<v.first<<" "<<neighbor.first<<" "<<std::endl;
        //          e_count ++;
        //      }
        // }
        // unsigned long long u_all=0,e_all=0;
        // MPI_Allreduce(&u_count, &u_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        // MPI_Allreduce(&e_count, &e_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        // int mpi_rank = havoqgt::comm_world().rank();
        // if(mpi_rank==0){
        //     std::cout<<"pruned graph:"<<u_all<<":"<<e_all<<std::endl;
        // }

        // open the result files
#ifdef ENABLE_CACHE
        auto alg_data = std::forward_as_tuple(pattern_graph, ordering, vertex_metadata, vertex_state_map, template_vertices, 
        plan, vertex_active_state_map, vertex_active_edges_map, result, offset_map, graph, 
        result_count, active_edge_map,time_col, stop, localCompute, query_cache);
#else
        auto alg_data = std::forward_as_tuple(pattern_graph, ordering, vertex_metadata, vertex_state_map, template_vertices, 
        plan, vertex_active_state_map, vertex_active_edges_map, result, offset_map, graph, 
        result_count, active_edge_map,time_col, stop, localCompute);
#endif
        
        typedef Enumeration_visitor<DataGraph> visitor_type;
        auto vq = create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(graph, alg_data);
        // auto vq = create_visitor_queue<visitor_type, Enumeration_queue_tds>(graph, alg_data);
        MPI_Barrier(MPI_COMM_WORLD);
        register_timer(&stop);

        std::vector<vertex_locator> start_visitors;
        // find the label of the first query vertex
        VertexData first_label = pattern_graph->vertex_data[ordering[0]];
        for (auto it = graph->vertices_begin(); it != graph->vertices_end(); it++){
            if (vertex_metadata[*it] == first_label){
                start_visitors.push_back(*it);
                // std::cout<<graph->locator_to_label(*it)<<" ";
            }
        }
        vq.init_visitor_traversal(start_visitors);
    }

    template <typename Key, typename Value>
    void gather_map(std::unordered_map<Key, Value>& mapping, int mpi_rank){
        int size = mapping.size();
        vector<std::pair<int, int>> buf_send, buf_rec;
        for(auto p : mapping){
            buf_send.push_back({p.first, p.second});
        }
        // if(mpi_rank == 0){
            // std::cout<<"send:"<<;
            // for(auto p : buf_send){
            //     std::cout<<"s:"<<mpi_rank<<":"<<p.first<<" "<<p.second<<std::endl;
            // }
        // }
        mpi_all_gather(buf_send, buf_rec, MPI_COMM_WORLD);
        mapping.clear();
        for(auto p : buf_rec){
            // if(mpi_rank == 1)
            //     std::cout<<mpi_rank<<" "<<p.first<<" "<<p.second<<std::endl;
            if(mapping.find(p.first) == mapping.end()){
                mapping.insert({p.first, p.second});
            }else{
                mapping[p.first] += p.second;
            }
        }
        // for(auto p: buf_send){
        //     if(mpi_rank == 0)
        //         std::cout<<mpi_rank<<" "<<p.first<<" "<<p.second<<std::endl;
        //     if(mapping.find(p.first) == mapping.end()){
        //         mapping.insert({p.first, p.second});
        //     }else{
        //         mapping[p.first] += p.second;
        //     }
        // }
        // if(mpi_rank == 0){
        //     for(auto p:mapping){
        //         std::cout<<"m:"<<p.first<<" "<<p.second<<std::endl;
        //     }
        // }
    }
    /**
 * VF3 ordering policy  
 */
    // generate the searching orderings
    template <typename PatternGraph, typename DataGraph, typename VertexMetadata, typename MPI_P>
    void generate_ordering_vf3(PatternGraph &pattern_graph, DataGraph *graph, VertexMetadata &vertex_metadata, MPI_P mpi_rank, MPI_P mpi_size, vector<Vertex> &ordering)
    {
        std::unordered_map<Vertex, double> Pf;
        std::unordered_map<VertexData, VertexData> classMap; //[label -> class]
        std::unordered_map<Vertex, Vertex> parentMap;        //[vertex -> parent of the vertex]
        std::unordered_map<int, std::unordered_map<VertexData, std::set<Vertex>>> P1set;
        // because of the graph being undirected, $S1set is depricated
        std::unordered_map<int, std::unordered_map<VertexData, std::set<Vertex>>> S1set;
        // std::vector<Vertex> NodeSequence;
        computeProbabilities(pattern_graph, graph, Pf, mpi_rank, mpi_size, vertex_metadata, classMap);
        // if(mpi_rank==1){
        //     std::cout<<"label frequency:"<<std::endl;
        //     for(auto v:Pf){
        //         std::cout<<v.first<<"->"<<v.second<<" ";
        //     }
        //     std::cout<<std::endl;
        // }
        generateNodeSequence(Pf, pattern_graph, ordering);
        // mpi_bcast(ordering, 0, MPI_COMM_WORLD);
    }

    // inDegreeMap is deprecated because of the graph being undirected
    template <typename DataGraph, typename Map, typename LabelMap, typename MPIRank,
              typename VertexMetadata>
    Count calculateDegreeMap(DataGraph *dataGraph, Map &inDegreeMap,
                             Map &outDegreeMap, LabelMap &labelMap,
                             MPIRank mpi_rank, MPIRank mpi_size,
                             VertexMetadata &vertex_metadata)
    {
        Count total_vertex_count = 0;
        for (auto vitr = dataGraph->vertices_begin(); vitr != dataGraph->vertices_end(); vitr++){
            // process the degree
            if ((uint32_t)mpi_rank ==(*vitr).owner()){ 
                // if the vertex is in the process
                Count out_de = dataGraph->outgoing_degree(*vitr);
                Count in_de = dataGraph->incoming_degree(*vitr);
                if (outDegreeMap.find(out_de) == outDegreeMap.end())
                {
                    outDegreeMap[out_de] = 1;
                }
                else
                {
                    outDegreeMap[out_de]++;
                }

                if (inDegreeMap.find(in_de) == inDegreeMap.end())
                {
                    inDegreeMap[in_de] = 1;
                }
                else
                {
                    inDegreeMap[in_de]++;
                }
            }
            // process the label
            auto label = vertex_metadata[*vitr];
            if (labelMap.find(label) == labelMap.end())
            {
                labelMap[label] = 1;
            }
            else
            {
                labelMap[label]++;
            }
        }
        gather_map(outDegreeMap, mpi_rank);
        gather_map(inDegreeMap, mpi_rank);
        gather_map(labelMap, mpi_rank);
        for(auto v:outDegreeMap){
            total_vertex_count += v.second;
        }
        // total_vertex_count = mpi_all_reduce(total_vertex_count, MPI_SUM, MPI_COMM_WORLD);
        // process the unordered_map into a struct list, in order to use the MPI send
        // outgoing degree, num
        // std::vector<std::pair<Count, Count>> sendBufOutDegree;
        // for (auto it = outDegreeMap.begin(); it != outDegreeMap.end(); it++)
        // {
        //     sendBufOutDegree.push_back(std::pair<Count, Count>(it->first, it->second));
        // }
        // // start to send the partial map
        // std::vector<std::pair<Count, Count>> recvBufOutDegree;
        // mpi_all_gather(sendBufOutDegree, recvBufOutDegree, MPI_COMM_WORLD); // every process will have the reduced results
        // // agregate the received results
        // outDegreeMap.clear();
        // for (auto it = recvBufOutDegree.begin(); it != recvBufOutDegree.end(); it++)
        // {
        //     if (outDegreeMap.find(it->first) == outDegreeMap.end())
        //     {
        //         outDegreeMap[it->first] = it->second;
        //     }
        //     else
        //     {
        //         outDegreeMap[it->first] += it->second;
        //     }
        //     total_vertex_count += it->second;
        // }

        // process the unordered_map into a struct list, in order to use the MPI send
        // incoming degree, num
        // std::vector<std::pair<Count, Count>> sendBufInDegree;
        // for (auto it = inDegreeMap.begin(); it != inDegreeMap.end(); it++)
        // {
        //     sendBufInDegree.push_back(std::pair<Count, Count>(it->first, it->second));
        // }
        // // start to send the partial map
        // std::vector<std::pair<Count, Count>> recvBufInDegree;
        // mpi_all_gather(
        //     sendBufInDegree, recvBufInDegree,
        //     MPI_COMM_WORLD); // every process will have the reduced results
        // // agregate the received results
        // inDegreeMap.clear();
        // for (auto it = recvBufInDegree.begin(); it != recvBufInDegree.end(); it++)
        // {
        //     if (inDegreeMap.find(it->first) == inDegreeMap.end())
        //     {
        //         inDegreeMap[it->first] = it->second;
        //     }
        //     else
        //     {
        //         inDegreeMap[it->first] += it->second;
        //     }
        // }

        // same procedure to process the label data
        // std::vector<std::pair<Count, Count>>
        //     sendBufLabel; // label(this type is supposed to $VertexDdata), num
        // for (auto it = labelMap.begin(); it != labelMap.end(); it++)
        // {
        //     sendBufLabel.push_back(std::pair<VertexData, Count>(it->first, it->second));
        // }
        // std::vector<std::pair<Count, Count>> recvBufLabel;
        // mpi_all_gather(sendBufLabel, recvBufLabel, MPI_COMM_WORLD);
        // labelMap.clear();
        // for (auto it = recvBufLabel.begin(); it != recvBufLabel.end(); it++)
        // {
        //     if (labelMap.find(it->first) == labelMap.end())
        //     {
        //         labelMap[it->first] = it->second;
        //     }
        //     else
        //     {
        //         labelMap[it->first] += it->second;
        //     }
        // }
        return total_vertex_count;
    }

    template <typename PatternGraph, typename DataGraph, typename PfMap,
              typename MPI_P, typename VertexMetadata, typename ClassMap>
    void computeProbabilities(PatternGraph &patternGraph, DataGraph *dataGraph,
                              PfMap &Pf, MPI_P mpi_rank, MPI_P mpi_size,
                              VertexMetadata &vertex_metadata, ClassMap &classMap)
    {
        std::unordered_map<Count, Count>
            out_degree_map; //[degree -> number of vertices]
        std::unordered_map<Count, Count>
            in_degree_map; //[degree -> number of vertices]
        std::unordered_map<VertexData, Count>
            label_map; //[VertexData -> number of vertices]
        // calculate all the maps
        Count total_vertex_count; // number of vertices in the data graph
        total_vertex_count =
            calculateDegreeMap(dataGraph, in_degree_map, out_degree_map, label_map,
                               mpi_rank, mpi_size, vertex_metadata);
        // if(mpi_rank==1){
        //     std::cout<<"degree:"<<std::endl;
        //     for(auto v:out_degree_map){
        //         std::cout<<v.first<<"->"<<v.second<<" ";
        //     }
        //     std::cout<<std::endl;
        //     std::cout<<"label:"<<std::endl;
        //     for(auto v:label_map){
        //         std::cout<<v.first<<"->"<<v.second<<" ";
        //     }
        //     std::cout<<std::endl;
        // }
        // iterate the pattern graph
        for (auto v = 0; v < patternGraph.vertex_data.size(); v++)
        {
            // process label
            double Pl = (double)label_map[patternGraph.vertex_data[v]] /
                        (double)total_vertex_count;
            // process degree
            Count patternDegree = patternGraph.vertex_degree[v];
            Count tmp = 0;

            for (auto it = out_degree_map.begin(); it != out_degree_map.end(); it++)
            {
                if (it->first >= patternDegree)
                {
                    tmp += it->second;
                }
            }
            double Pd = (double)tmp / (double)total_vertex_count;
            Pf[v] = Pl * Pd;
        }

        // process the classMap, by default each label maps to a unique class
        for (auto it = label_map.begin(); it != label_map.end(); it++)
        {
            classMap[it->first] = it->first;
        }
    }

    struct BindData
    {
        Vertex v;
        int dM; // Node mapping degree: the number of incoming and outgoing edges
                // between u and all the nodes that are already inside NG1;
        int degree;
        double probability;
        BindData(Vertex vertex, int dm, int de, double pr)
            : v(vertex), dM(dm), degree(de), probability(pr) {}
    };

    /* Generate ordering by first comparing dM, then Pf, finally degree of the nodes
 * in a pattern graph. */
    template <typename PfMap, typename PatternGraph, typename NodeVector>
    void generateNodeSequence(PfMap &Pf, PatternGraph &patternGraph,
                              NodeVector &result)
    {
        std::vector<BindData> tmp;
        for (auto v = 0; v < patternGraph.vertex_degree.size(); v++)
        {
            tmp.push_back(BindData(v, 0, patternGraph.vertex_degree[v], Pf[v]));
        }
        for (int j = 0; j < patternGraph.vertex_degree.size(); j++)
        {
            BindData current = tmp[0];
            auto currentIndex = tmp.begin();

            // Select a node of the pattern graph from tmp, and insert into result
            for (auto it = tmp.begin(); it != tmp.end(); it++)
            {
                BindData t = *it;
                if (current.dM < t.dM)
                {
                    current = t;
                    currentIndex = it;
                }
                else if (current.dM == t.dM)
                {
                    if (current.probability > t.probability)
                    {
                        current = t;
                        currentIndex = it;
                    }
                    else if (current.probability == t.probability)
                    {
                        if (current.degree < t.degree)
                        {
                            current = t;
                            currentIndex = it;
                        }
                    }
                }
            }
            result.push_back(current.v);
            tmp.erase(currentIndex);

            // Update the value of dM for the node in tmp
            for (auto it = tmp.begin(); it != tmp.end(); it++)
            {
                for (auto it1 = result.begin(); it1 != result.end(); it1++)
                {
                    for (auto e = patternGraph.vertices[it->v];
                         e < patternGraph.vertices[it->v + 1]; e++)
                    {
                        if (patternGraph.edges[e] ==
                            *it1)
                        { // vertex it and it1 are neighbours
                            it->dM++;
                            break;
                        }
                    }
                }
            }
        }
    }

    // because of the graph being undirected, thus Set S is depricated
    template <typename PatternGraph, typename NodeSequence, typename ClassMap,
              typename ParentMap, typename PS>
    void classifyNodes(PatternGraph &patternGraph, NodeSequence &nodeList,
                       ClassMap &classMap, ParentMap &parentMap, PS &P, PS &S)
    {
        for (int level = 1; level <= patternGraph.vertex_count; level++)
        {
            std::unordered_map<VertexClass, std::set<Vertex>> Pset;
            std::unordered_map<VertexClass, std::set<Vertex>> Sset;
            for (auto it = nodeList.begin();
                 std::distance(nodeList.begin(), it) < level; it++)
            {
                Vertex v = *it;
                for (auto e = patternGraph.vertices[v]; e < patternGraph.vertices[v + 1];
                     e++)
                { // neighbors of the v
                    Vertex neighbor = patternGraph.edges[e];
                    VertexClass classes = classMap[patternGraph.vertex_data[neighbor]];
                    auto iter = std::find(nodeList.begin(), nodeList.end(), neighbor);
                    if (std::distance(nodeList.begin(), iter) >= level)
                    {
                        if (Pset.find(classes) == Pset.end())
                        {
                            std::set<Vertex> tmp;
                            Pset[classes] = tmp;
                        }
                        Pset[classes].insert(neighbor);
                        if (parentMap.find(neighbor) == parentMap.end())
                        {
                            parentMap[neighbor] = nodeList[level - 1];
                        }
                    }
                }
            }
            P[level] = Pset;
        }
    }
} // namespace prunejuice
