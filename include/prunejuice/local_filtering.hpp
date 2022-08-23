#pragma once

#include <bitset>
#include <boost/dynamic_bitset.hpp>

#include <havoqgt/detail/visitor_priority_queue.hpp>
#include <havoqgt/visitor_queue.hpp>

#define MAX_PATTERN_SIZE 32

#define OUTPUT_RESULT

namespace prunejuice
{

  template <typename Visitor>
  class ldf_queue
  {
  public:
    ldf_queue() {}

    bool push(Visitor const &element)
    {
      data.push_back(element);
      return true;
    }

    void pop() { data.pop_back(); }

    Visitor const &top() { return data.back(); }

    size_t size() const { return data.size(); }

    bool empty() const { return data.empty(); }

    void clear() { data.clear(); }

  protected:
    std::vector<Visitor> data;
  }; // end of class

  template <typename Graph, typename Vertex, typename VertexData, typename BitSet,
            typename NLFMap>
  class ldf_visitor
  {
    typedef typename Graph::vertex_locator vertex_locator;
    typedef typename Graph::edge_iterator eitr_type;
    typedef prunejuice::vertex_state<Vertex, VertexData, BitSet, NLFMap> VertexState;

  public:
    ldf_visitor() : msg_type(0) {}

    ldf_visitor(vertex_locator _vertex) : vertex(_vertex), msg_type(0) {}

    ldf_visitor(vertex_locator _vertex, vertex_locator parent_, VertexData parent_label_, uint8_t msg_type_)
        : vertex(_vertex), parent(parent_), parent_label(parent_label_), msg_type(msg_type_) {}

    // ldf_visitor(vertex_locator _vertex, vertex_locator _parent,
    //             uint8_t _msg_type = 0)
    //     : vertex(_vertex), parent(_parent), msg_type(_msg_type) {}

    // ldf_visitor(vertex_locator _vertex, vertex_locator _parent,
    //             BitSet _parent_template_neighbors, uint8_t _msg_type)
    //     : vertex(_vertex),
    //       parent(_parent),
    //       parent_template_neighbors(_parent_template_neighbors),
    //       msg_type(_msg_type) {}

    // ldf_visitor(vertex_locator _vertex, vertex_locator _parent,
    //             BitSet _parent_template_neighbors, VertexData _parent_label,
    //             uint8_t _msg_type)
    //     : vertex(_vertex),
    //       parent(_parent),
    //       parent_template_neighbors(_parent_template_neighbors),
    //       parent_label(_parent_label),
    //       msg_type(_msg_type) {}

    ~ldf_visitor() {}

    template <typename AlgData>
    bool pre_visit(AlgData &alg_data) const{
      return true;
    }

    template <typename VisitorQueueHandle, typename AlgData>
    bool init_visit(Graph &g, VisitorQueueHandle vis_queue,
                    AlgData &alg_data) const
    {
      return visit(g, vis_queue, alg_data);
    }

    template <typename VisitorQueueHandle, typename AlgData>
    bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data) const{
      auto & global_init_step = std::get<7>(alg_data);
      auto & super_var = std::get<8>(alg_data);
      auto & vertex_metadata = std::get<1>(alg_data);
      auto & vertex_state_map = std::get<3>(alg_data);
      auto & active_vertex = std::get<4>(alg_data);
      auto & pattern_graph = std::get<2>(alg_data);
      auto & update_active_edges = std::get<9>(alg_data);
      auto mpi_rank = havoqgt::comm_world().rank();
      if(update_active_edges){
        // std::cout<<"acc:"<<g.locator_to_label(vertex)<<std::endl;
        auto itf = vertex_state_map.find(g.locator_to_label(vertex));
        if(itf == vertex_state_map.end()){
          return false;
        }

        if(msg_type == 0){
          bool pruned = false;
          for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
            vertex_locator neighbor = eitr.target();
            //std::cout << "Visiting neighbor: " << g.locator_to_label(neighbor) << std::endl;
            ldf_visitor new_visitor(neighbor, vertex, vertex_metadata[vertex], 1);
            // std::cout<<"send:"<<g.locator_to_label(vertex)<<"->"<<g.locator_to_label(neighbor)<<std::endl;
            vis_queue->queue_visitor(new_visitor);
          }
        }else{
          // update the active_edges
          auto& vertex_active_edges_map = std::get<5>(alg_data);
          // std::cout<<"rec:"<<g.locator_to_label(vertex)<<"<-"<<g.locator_to_label(parent)<<std::endl;
          auto itf = vertex_active_edges_map[vertex].find(g.locator_to_label(parent));
          bool r = (itf != vertex_active_edges_map[vertex].end());
          auto insert_status = vertex_active_edges_map[vertex].insert({g.locator_to_label(parent), 1});
          // if(!insert_status.second){
          //   if(r){
          //     std::cout<<"ERROR: already exists"<<std::endl;
          //   }else{
          //     std::cout<<"ERROR: updating the active edges"<<std::endl;
          //   }
          // }
        }
        return true;
      }

      if(global_init_step && super_var == 0){
        std::get<10>(alg_data) = 0;
        if(msg_type == 0){
          for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
            vertex_locator neighbor = eitr.target();
            // if(g.locator_to_label(neighbor) == 2206366){
            //     std::cout<<"send:"<<mpi_rank<<":"<<vertex.owner()<<":"<<g.locator_to_label(vertex)<<"->"<<g.locator_to_label(neighbor)<<":"<<neighbor.owner()<<std::endl;
            //     // for(auto ii=g.edges_begin(vertex); ii!=g.edges_end(vertex); ++ii){
            //     //    std::cout<<g.locator_to_label(ii.target())<<",";
            //     // }
            //     // std::cout<<"}"<<std::endl;
            // }
            //std::cout << "Visiting neighbor: " << g.locator_to_label(neighbor) << std::endl;
            ldf_visitor new_visitor(neighbor, vertex, vertex_metadata[vertex], 1);
            vis_queue->queue_visitor(new_visitor);
          }
        }else{
          // initialise the basic data structures
          // if(g.locator_to_label(vertex) == 2206366){
          //     std::cout<<"recv:"<<mpi_rank<<":"<<vertex.owner()<<":"<<g.locator_to_label(vertex)<<"<-"<<g.locator_to_label(parent)<<":"<<parent.owner()<<std::endl;
          // }
          auto itf = vertex_state_map.find(g.locator_to_label(vertex));
          if(itf == vertex_state_map.end()){
            auto insert_status = vertex_state_map.insert({g.locator_to_label(vertex), VertexState()});
            if(!insert_status.second){
              std::cout<<"Error updating the vertex_state_map"<<std::endl;
            }
            itf = insert_status.first;
          }
          // candidate set
          itf->second.template_vertices.set(); // unknown property
          itf->second.template_neighbors.set();
          // construct the NLF
          auto itf1 = itf->second.neighbor_label_frequency_map.find(parent_label);
          if(itf1 == itf->second.neighbor_label_frequency_map.end()){
            itf->second.neighbor_label_frequency_map.insert({parent_label, 1});
          }else{
            itf1->second ++;
          }
        }
      }else{
        // start pruning process
        auto itf = vertex_state_map.find(g.locator_to_label(vertex));
        if(itf == vertex_state_map.end()){
          return false;
        }

        if(msg_type == 0){
          // checking the LDF
          // std::cout<<"test:"<<g.locator_to_label(vertex)<<":t"<<std::endl;

          auto label = vertex_metadata[vertex];
          for(int i=0;i<itf->second.template_neighbors.size(); ++i){
            if(i >= pattern_graph.vertex_data.size()){
              itf->second.template_neighbors.reset(i);
            }else if(label != pattern_graph.vertex_data[i]){
              itf->second.template_neighbors.reset(i);
            }
          }
          // if(label!=13&&label!=15&&label!=12&&label!=9&&label!=16&&!itf->second.template_neighbors.none()){
          //     std::cout<<g.locator_to_label(vertex)<<":"<<label<<":"<<itf->second.template_neighbors<<std::endl;
          //     exit(0);
          // }
          // std::cout<<"test:"<<g.locator_to_label(vertex)<<":s:"<<itf->second.template_neighbors.size()<<std::endl;
          // checking NLF
          for(int i=0;i<pattern_graph.vertex_data.size(); ++i){
            if(itf->second.template_neighbors.test(i)){
              // std::cout<<"cc:"<<g.locator_to_label(vertex)<<":"<<i<<":"<<itf->second.template_neighbors.size()<<std::endl;
              for(auto &n : pattern_graph.vertex_neighbor_data_count_map[i]){
                auto itf1 = itf->second.neighbor_label_frequency_map.find(n.first);
                if(itf1 == itf->second.neighbor_label_frequency_map.end()){
                  itf->second.template_neighbors.reset(i);
                }else if(itf1->second < n.second){
                  itf->second.template_neighbors.reset(i);
                }
              }
              // std::cout<<"cc1:"<<g.locator_to_label(vertex)<<":"<<i<<std::endl;
            }
          }
          // if(label!=13&&label!=15&&label!=12&&label!=9&&label!=16&&!itf->second.template_neighbors.none()){
          //     std::cout<<g.locator_to_label(vertex)<<":a:"<<label<<":"<<itf->second.template_neighbors<<std::endl;
          //     exit(0);
          // }
          // std::cout<<"test:"<<g.locator_to_label(vertex)<<":o"<<std::endl;
          if(itf->second.template_neighbors.none()){
            std::get<10>(alg_data) = 0;
            // remove the vertex
            vertex_state_map.erase(itf);
            // if(vertex_state_map.find(g.locator_to_label(vertex)) != vertex_state_map.end()){
            //     std::cout<<g.locator_to_label(vertex)<<":error:"<<label<<":"<<itf->second.template_neighbors<<std::endl;
            //     exit(0);
            // }
            // broadcast the updates
            for(auto eitr = g.edges_begin(vertex); eitr != g.edges_end(vertex); ++eitr) {
              vertex_locator neighbor = eitr.target();
              // if(g.locator_to_label(neighbor) == 2206366 || g.locator_to_label(vertex)==2206366){
              //   std::cout<<"send:"<<mpi_rank<<":"<<vertex.owner()<<":"<<g.locator_to_label(vertex)<<"->"<<g.locator_to_label(neighbor)<<":"<<neighbor.owner()<<std::endl;
              // }
              // std::cout<<"pruned:"<<g.locator_to_label(vertex)<<"->"<<g.locator_to_label(neighbor)<<std::endl;
              //std::cout << "Visiting neighbor: " << g.locator_to_label(neighbor) << std::endl;
              ldf_visitor new_visitor(neighbor, vertex, vertex_metadata[vertex], 1);
              vis_queue->queue_visitor(new_visitor);
            }
          }
        }else{
          // // update the NLF
          //  if(g.locator_to_label(vertex) == 2206366 || g.locator_to_label(parent)==2206366){
          //      std::cout<<"recv:"<<mpi_rank<<":"<<vertex.owner()<<":"<<g.locator_to_label(vertex)<<"<-"<<g.locator_to_label(parent)<<":"<<parent.owner()<<std::endl;
          //  }
          auto itf = vertex_state_map.find(g.locator_to_label(vertex));
          // std::cout<<"update:"<<g.locator_to_label(vertex)<<"<-"<<g.locator_to_label(parent)<<std::endl;
          auto itf1 = itf->second.neighbor_label_frequency_map.find(parent_label);
          if(itf1 == itf->second.neighbor_label_frequency_map.end()){
            std::cout<<"ERROR: the entry should exists"<<std::endl;
          }else{
            itf1->second --;
          }
        }
      }
      return true;
    }

    friend inline bool operator>(const ldf_visitor &v1, const ldf_visitor &v2)
    {
      return v1.msg_type < v2.msg_type;
    }

    vertex_locator vertex;
    vertex_locator parent; // parent is a neighbor
    VertexData parent_label;
    uint8_t msg_type;
  }; // end of class




// template <typename Key, typename Value>
// void gather_map(std::unordered_map<Key, Value>& mapping, int mpi_rank){
//         int size = mapping.size();
//         vector<std::pair<int, int>> buf_send, buf_rec;
//         for(auto p : mapping){
//             buf_send.push_back({p.first, p.second});
//         }
//         mpi_all_gather(buf_send, buf_rec, MPI_COMM_WORLD);
//         mapping.clear();
//         for(auto p : buf_rec){
//             if(mapping.find(p.first) == mapping.end()){
//                 mapping.insert({p.first, p.second});
//             }else{
//                 mapping[p.first] += p.second;
//             }
//         }
// }     

template <typename Vertex, typename VertexData, typename VertexMetaData,
            typename BitSet, typename NLFMap, typename TGraph,
            typename PatternGraph, typename VertexStateMap, typename VertexActive,
            typename VertexUint8MapCollection, typename TemplateVertex>
  void label_propagation_label_degree_filtering_bsp(
      TGraph *g, VertexMetaData &vertex_metadata, PatternGraph &pattern_graph,
      VertexStateMap &vertex_state_map, VertexActive &vertex_active,
      VertexUint8MapCollection &vertex_active_edges_map,
      TemplateVertex &template_neighbors, bool global_init_step,
      bool &global_not_finished, std::ofstream &superstep_result_file,
      std::ofstream &message_count_result_file,
      std::ofstream &active_vertices_count_result_file,
      std::ofstream &active_edges_count_result_file)
  {
    int mpi_rank = havoqgt::comm_world().rank();
    uint64_t superstep_var = 0;
    uint64_t &superstep_ref = superstep_var;
    bool update_active_edges = false;
    int not_finished = false;

    typedef ldf_visitor<TGraph, Vertex, VertexData, BitSet, NLFMap> visitor_type;
    auto alg_data = std::forward_as_tuple(
        g, vertex_metadata, pattern_graph, vertex_state_map, vertex_active,
        vertex_active_edges_map, template_neighbors, global_init_step,
        superstep_var, update_active_edges, not_finished);
    auto vq =
        havoqgt::create_visitor_queue<visitor_type, havoqgt::detail::visitor_priority_queue>(g, alg_data);

    if (mpi_rank == 0)
    {
      std::cout << "Local Constraint Checking ... " << std::endl;
    }
    if(global_init_step){
        for(auto v:vertex_state_map){
          vertex_active_edges_map[g->label_to_locator(v.first)].clear();
          std::cout<<v.first<<":{";
          for(auto l : vertex_active_edges_map[g->label_to_locator(v.first)]){
            std::cout<<l.first<<",";
          }
          std::cout<<"}"<<std::endl;
        }
    }
    // init
    if(global_init_step){
        double time_start = MPI_Wtime();
        vq.init_visitor_traversal();
        MPI_Barrier(MPI_COMM_WORLD);
        double time_end = MPI_Wtime();
        if(mpi_rank==0)
            std::cout<<"init time:"<<time_end-time_start<<std::endl;
        ++ superstep_var;    
    }else{
        unsigned long long after_size=0, total_size;
        for(auto&v : vertex_state_map){
            ++after_size;   
        }
      MPI_Allreduce(&after_size, &total_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      if(mpi_rank==0){
          std::cout<<"before"<<superstep_var-1<<":"<<total_size<<std::endl;
      }

    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(int i=0;i<10 && not_finished != true;++i){
      // initialise start vertices
      std::vector<vertex_locator> start_visitors;
      for(auto& v: vertex_state_map){
          start_visitors.push_back(g->label_to_locator(v.first));
      }
      not_finished = 1;
      double time_start = MPI_Wtime();
      vq.init_visitor_traversal();
      MPI_Barrier(MPI_COMM_WORLD);
      double time_end = MPI_Wtime();
      int local_not_finished = not_finished;
      MPI_Allreduce(&local_not_finished, &not_finished, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      // not_finished = havoqgt::mpi_all_reduce(not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);

      if (mpi_rank == 0){
        std::cout << " | Time : " << time_end - time_start << std::endl;
      }

      if (mpi_rank == 0){
        std::cout << "Local Constraint Checking | Local Finished Status : ";
        if (not_finished){
          std::cout << "Continue" << std::endl;
        }
        else{
          std::cout << "Stop" << std::endl;
        }
      }
      ++ superstep_var;
      for(auto vitr = g->vertices_begin(); vitr != g->vertices_end(); ++vitr){
          if(vertex_state_map.find(g->locator_to_label(*vitr)) != vertex_state_map.end()){
              vertex_active[*vitr] = true;
          }else{
              vertex_active[*vitr] = false;
          }
      } 
      unsigned long long after_size=0, total_size;
      for(auto&v : vertex_state_map){
          // initialise the vertex_active
          // if(loc.is_delegate()){
          //     std::cout<<"has delegete:"<<mpi_rank<<":"<<v.first<<std::endl;
          // }
          // if(label!=13&&label!=15&&label!=12&&label!=9&&label!=16&&i>0){
          //     auto itf = vertex_state_map.find(v.first);
          //     std::cout<<v.first<<":"<<mpi_rank<<":"<<loc.is_delegate()<<":"<<g->master(loc)<<":"<<":al:"<<label<<":"<<itf->second.template_neighbors<<std::endl;
          //     exit(0);
          // }
          ++after_size;   
      }
      MPI_Allreduce(&after_size, &total_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      if(mpi_rank==0){
          std::cout<<"iteration:"<<superstep_var-1<<":"<<total_size<<std::endl;
      }
    }
    global_init_step = false;
    update_active_edges = true;
    
    vq.init_visitor_traversal();

    // updating the neighbor_candidates
    for(auto v:vertex_state_map){
      template_neighbors[g->label_to_locator(v.first)] = v.second.template_neighbors.to_ulong();
    }

    // for(auto v:vertex_state_map){
    //   if(v.first == 0){
    //     std::cout<<mpi_rank<<":"<<v.first<<":{";
    //     for(auto l : v.second.neighbor_label_frequency_map){
    //       std::cout<<l.first<<":"<<l.second<<",";
    //     }
    //     std::cout<<"}"<<std::endl;
    //   }
    // }
    // for(auto v:vertex_state_map){
    //   // auto itf = vertex_active_edges_map.find(g->label_to_locator(v.first));
    //   std::cout<<v.first<<":{";
    //   for(auto l : vertex_active_edges_map[g->label_to_locator(v.first)]){
    //     std::cout<<l.first<<",";
    //   }
    //   std::cout<<"}"<<std::endl;
    // }
    return;
  }
} // namespace prunejuice
