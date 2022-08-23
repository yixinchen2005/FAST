// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#include <bitset>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <set>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/interprocess/managed_heap_memory.hpp>

#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
///#include <havoqgt/environment.hpp>

#include <metadata/vertex_data_db.hpp>
#include <metadata/vertex_data_db_degree.hpp>

// #include <prunejuice/template.hpp>
#include <prunejuice/config.hpp>
#include <prunejuice/template_refine.hpp>
#include <prunejuice/non_local_constraint.hpp>
#include <prunejuice/algorithm_state.hpp>
#include <prunejuice/local_constraint_checking.hpp>
#include <prunejuice/non_local_constraint_checking.hpp>
#include <prunejuice/non_local_constraint_checking_unique.hpp>
#include <prunejuice/non_local_constraint_checking_tds_batch.hpp> // TODO: optimize for batch_size = mpi_size

#include <prunejuice/enumerate.hpp>
#include <prunejuice/query_cache.hpp>
#include <prunejuice/random_sample.hpp>

#ifdef ENABLE_MY_LCC
#include <prunejuice/local_filtering.hpp>
#endif

//#define OUTPUT_RESULT
//#define ENABLE_BLOCK

#define TP_ASYNC
//#define TP_BATCH

///namespace hmpi = havoqgt::mpi;
///using namespace havoqgt::mpi;
using namespace havoqgt;
//using namespace prunejuice;

///typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
//typedef havoqgt::delegate_partitioned_graph<segment_manager_t> graph_type;
typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

template <typename T>
using DelegateGraphVertexDataSTDAllocator = graph_type::vertex_data<T, std::allocator<T>>;

template <typename T>
using DelegateGraphEdgeDataSTDAllocator = graph_type::edge_data<T, std::allocator<T>>;


template<typename Graph, typename VertexMap, typename EdgeMap>
void get_pruned_graph_size(Graph* graph, VertexMap& vertex_state_map, EdgeMap& vertex_active_edges_map, unsigned long long&remain_all_vertices, unsigned long long& remain_all_edges){
  MPI_Barrier(MPI_COMM_WORLD); 
  unsigned long long remain_vertices=0,remain_edges=0;
  for (auto &v : vertex_state_map)
  {
    auto v_locator = graph->label_to_locator(v.first);
    ++ remain_vertices;
    for (auto &n : vertex_active_edges_map[v_locator])
    {
      ++ remain_edges;
    }
  }
  MPI_Allreduce(&remain_vertices, &remain_all_vertices, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&remain_edges, &remain_all_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
}

void usage()
{
  //if(havoqgt_env()->world_comm().rank() == 0) {
  if (comm_world().rank() == 0)
  {
    std::cerr << "Usage: -i <string> -p <string> -o <string>\n"
              << " -i <string>   - input graph base filename (required)\n"
              << " -b <string>   - backup graph base filename. If set, \"input\" graph will be deleted if it exists\n"
              << " -v <string>   - vertex metadata base filename (optional, Default is degree based metadata)\n"
              << " -e <string>   - edge metadata base filename (optional)\n"
              << " -p <string>   - pattern base directory (required)\n"
              << " -o <string>   - output base directory (required)\n"
              //<< " -x <int>      - Token Passing batch size (optional, Default/max batch size is "
              //  << havoqgt_env()->world_comm().size() << ", Min batch size is 1)\n"
              //  << comm_world().size() << " , min batch size is 1)\n"
              << " -x <int>      - Token Passing batch count (optional, Default/min batch count is 1, max batch count is "
              << comm_world().size() << "\n"
              << " -h            - print help and exit\n\n";
  }
}



void parse_cmd_line(int argc, char **argv, std::string &graph_input,
                    std::string &backup_graph_input, std::string &vertex_metadata_input,
                    std::string &edge_metadata_input, std::string &pattern_input,
                    std::string &result_output, uint64_t &tp_vertex_batch_size)
{

  //if(havoqgt_env()->world_comm().rank() == 0) {
  if (comm_world().rank() == 0)
  {
    std::cout << "CMD Line :";
    for (int i = 0; i < argc; ++i)
    {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool print_help = false;
  std::bitset<3> required_input;
  required_input.reset();

  char c;
  while ((c = getopt(argc, argv, "i:b:v:e:p:o:x:h ")) != -1)
  {
    switch (c)
    {
    case 'h':
      print_help = true;
      break;
    case 'i':
      graph_input = optarg;
      required_input.set(0);
      break;
    case 'b':
      backup_graph_input = optarg;
      break;
    case 'v':
      vertex_metadata_input = optarg;
      break;
    case 'e':
      edge_metadata_input = optarg;
      break;
    case 'p':
      pattern_input = optarg;
      required_input.set(1);
      break;
    case 'o':
      result_output = optarg;
      required_input.set(2);
      break;
    case 'x':
      tp_vertex_batch_size = std::stoull(optarg);
      if (tp_vertex_batch_size < 1 || tp_vertex_batch_size > comm_world().size())
      {
        print_help = true;
      }
      else if (tp_vertex_batch_size > 1)
      {
        tp_vertex_batch_size = comm_world().size() / tp_vertex_batch_size;
      }
      else
      {
        tp_vertex_batch_size = comm_world().size();
      }
      break;
    default:
      std::cerr << "Unrecognized Option : " << c << ", Ignore." << std::endl;
      print_help = true;
      break;
    }
  }

  if (print_help || !required_input.all())
  {
    usage();
    exit(-1);
  }
}

int main(int argc, char **argv)
{
  //typedef hmpi::delegate_partitioned_graph<segment_manager_t> graph_type;
  //typedef hmpi::delegate_partitioned_graph
  //  <typename segment_manager_t::template allocator<void>::type> graph_type;

  int mpi_rank(0), mpi_size(0);

  // havoqgt_init
  //havoqgt::havoqgt_init(&argc, &argv);
  havoqgt::init(&argc, &argv);
  {

    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    ///havoqgt::get_environment();

    if (mpi_rank == 0)
    {
      std::cout << "MPI Initialized With " << mpi_size << " Ranks." << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_graph_input;
    std::string vertex_metadata_input;
    std::string edge_metadata_input;
    std::string pattern_input;
    std::string result_output;
    ///uint64_t tp_vertex_batch_size = havoqgt_env()->world_comm().size();
    uint64_t tp_vertex_batch_size = comm_world().size();

    parse_cmd_line(argc, argv, graph_input, backup_graph_input,
                   vertex_metadata_input, edge_metadata_input, pattern_input, result_output,
                   tp_vertex_batch_size);

    std::string pattern_dir = pattern_input;
    std::string result_dir = result_output;

    MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////////////////////////////////////////////////////////////

    // load graph

    if (mpi_rank == 0)
    {
      std::cout << "Loading Graph ... " << std::endl;
    }

    if (backup_graph_input.size() > 0)
    {
      distributed_db::transfer(backup_graph_input.c_str(), graph_input.c_str());
    }

    distributed_db ddb(db_open_read_only(), graph_input.c_str());

    auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    typedef uint8_t edge_data_type;
    typedef graph_type::edge_data<edge_data_type, distributed_db::allocator<edge_data_type>> edge_data_t;

    auto edge_data_ptr = ddb.get_manager()->find<edge_data_t>("graph_edge_data_obj").first;
    //  assert(edge_data_ptr != nullptr);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0)
    {
      std::cout << "Done Loading Graph." << std::endl;
    }

    /////////////////////////////////////////////////////////////////////////////

    // pattern matching
    {

      // types used by the delegate partitioned graph
      typedef typename graph_type::vertex_iterator vitr_type;
      typedef typename graph_type::vertex_locator vloc_type;
      //typedef typename graph_type::edge_iterator eitr_type;

      typedef uint64_t Vertex;
      typedef uint64_t Edge;
      typedef uint64_t VertexData; 
      typedef uint8_t EdgeData;

      typedef uint64_t VertexRankType; // TODO: delete

      static constexpr size_t max_bit_vector_size = 32; // TODO:
      static constexpr size_t max_template_vertex_count = 32;
      typedef std::bitset<max_bit_vector_size> BitSet; // TODO: rename to TemplateVertexSet
      typedef BitSet TemplateVertexBitSet;
      typedef uint16_t TemplateVertexType; // TODO: rename to TemplateVertexBitSetToUint

      typedef uint8_t Boolean; // TODO: replace all bool with Boolean?

      // TODO: mmap
      //typedef graph_type::vertex_data<VertexData, SegmentAllocator<VertexData> > VertexMetadata;
      //typedef graph_type::vertex_data<VertexRankType, SegmentAllocator<VertexRankType> > VertexRank;
      //typedef graph_type::vertex_data<bool, SegmentAllocator<bool> > VertexActive;
      //typedef graph_type::vertex_data<uint64_t, SegmentAllocator<uint64_t> > VertexIteration;

      typedef graph_type::vertex_data<VertexData, std::allocator<VertexData>> VertexMetadata;
      typedef graph_type::vertex_data<Boolean, std::allocator<Boolean>> VertexActive;                         // TODO: solution_graph // TODO: you are mixing bool and uint!
      typedef graph_type::vertex_data<TemplateVertexType, std::allocator<TemplateVertexType>> TemplateVertex; // TODO: solution_graph, rename to VertexTemplateVertexBitSetToUint

      typedef graph_type::vertex_data<uint64_t, std::allocator<uint64_t>> VertexIteration;        // TODO: delete
      typedef graph_type::vertex_data<VertexRankType, std::allocator<VertexRankType>> VertexRank; // TODO: delete

      //typedef vertex_state<uint8_t> VertexState;
      ///  typedef prunejuice::vertex_state_generic<Vertex, VertexData, uint8_t, BitSet> VertexState;
      typedef std::unordered_map<VertexData, uint64_t> NLFMap;
#ifdef ENABLE_MY_LCC
      typedef prunejuice::vertex_state<Vertex, VertexData, BitSet, NLFMap> VertexState;
#else
      typedef prunejuice::vertex_state<Vertex, VertexData, BitSet> VertexState;
#endif
      typedef std::unordered_map<Vertex, VertexState> VertexStateMap; // TODO: solution_graph

      typedef std::unordered_set<Vertex> VertexSet;
      typedef graph_type::vertex_data<VertexSet, std::allocator<VertexSet>> VertexSetCollection;

      typedef std::unordered_map<Vertex, uint8_t> VertexUint8Map;
      typedef graph_type::vertex_data<VertexUint8Map, std::allocator<VertexUint8Map>> VertexUint8MapCollection;

      typedef graph_type::edge_data<EdgeData, std::allocator<EdgeData>> EdgeMetadata;
      typedef graph_type::edge_data<Boolean, std::allocator<Boolean>> EdgeActive; // TODO: solution_graph

      typedef std::vector<Boolean> VectorBoolean;

      //////////////////////////////////////////////////////////////////////////////

      if (mpi_rank == 0)
      {
        std::cout << "Pattern Matching ... " << std::endl;
      }

      //////////////////////////////////////////////////////////////////////////////

      double time_start = MPI_Wtime();
      double time_end = MPI_Wtime();

      // per rank containers
      VertexStateMap vertex_state_map;

      // vertex containers

      // TODO: need a new alloc_inst to use bip/mmap
      //  VertexMetadata vertex_metadata(*graph, alloc_inst);
      //  VertexRank vertex_rank(*graph, alloc_inst);
      //  VertexActive vertex_active(*graph, alloc_inst);
      //  VertexIteration vertex_iteration(*graph, alloc_inst);

      VertexMetadata vertex_metadata(*graph);
      VertexActive vertex_active(*graph);
      TemplateVertex template_vertices(*graph);
      VertexUint8MapCollection vertex_active_edges_map(*graph);
      VertexSetCollection vertex_token_source_set(*graph); // per vertex set
                                                     //--  VertexSetCollection token_source_edge_set(*graph); // per vertex set // edge aware
      // std::cout<<template_vertices[graph->label_to_locator(0)]<<std::endl;
      BitSet vertex_template_vertices(template_vertices[graph->label_to_locator(0)]);
      // std::cout<<vertex_active_edges_map<<std::endl;
      //  VertexRank vertex_rank(*graph);
      uint8_t vertex_rank;      // TODO: dummy
                                //  VertexIteration vertex_iteration(*graph);
      uint8_t vertex_iteration; // TODO: dummy

      // edge containers
      //EdgeMetadata edge_metadata(*graph);
      ////EdgeActive edge_active(*graph);

      if (mpi_rank == 0)
      {
        std::cout << "Pattern Matching | Allocated Vertex and Edge Containers"
                  << std::endl;
      }

      /////////////////////////////////////////////////////////////////////////////

      // application parameters // TODO: command line input

      // write vertex data to file
      bool do_output_vertex_data = true; // TODO: ?

      size_t token_passing_algo = 0; // TODO: ?

      MPI_Barrier(MPI_COMM_WORLD);

      /////////////////////////////////////////////////////////////////////////////

      // build the distributed vertex data db
      time_start = MPI_Wtime();

      //if (use_degree_as_vertex_data) {
      if (vertex_metadata_input.size() > 0)
      {
        vertex_data_db_nostdfs<graph_type, VertexMetadata, Vertex, VertexData>
            //(graph, vertex_metadata, vertex_data_input_filename, 10000);
            (graph, vertex_metadata, vertex_metadata_input, 10000);
        // TODO: each rank reads 10K lines from file at a time
      }
      else
      {
        vertex_data_db_degree<graph_type, VertexMetadata, Vertex, VertexData>(graph, vertex_metadata);
      }

      MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this?
      time_end = MPI_Wtime();
      if (mpi_rank == 0)
      {
        std::cout << "Pattern Matching Time | Vertex Data DB : "
                  << time_end - time_start << std::endl;
      }

      if (do_output_vertex_data)
      {
        std::string vertex_data_filename = result_dir + //"/" + std::to_string(0) +
                                           "/all_ranks_vertex_data/vertex_data_" + std::to_string(mpi_rank);
        std::ofstream vertex_data_file(vertex_data_filename, std::ofstream::out);

        for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();
             ++vitr)
        {
          vloc_type vertex = *vitr;
          vertex_data_file << mpi_rank << ", l, " << graph->locator_to_label(vertex)
                           << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
        }

        for (vitr_type vitr = graph->delegate_vertices_begin();
             vitr != graph->delegate_vertices_end(); ++vitr)
        {
          vloc_type vertex = *vitr;

          if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank))
          {
            vertex_data_file << mpi_rank << ", c, " << graph->locator_to_label(vertex)
                             << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
          }
          else
          {
            vertex_data_file << mpi_rank << ", d, " << graph->locator_to_label(vertex)
                             << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";
          }
        }

        vertex_data_file.close();
      }
      // collect edge and vertex count before pruning
      unsigned long long vertex_count=0, edge_count=0;
      for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end();++vitr){
        ++vertex_count;
        for(auto eitr = graph->edges_begin(*vitr); eitr != graph->edges_end(*vitr); ++eitr) {
          ++ edge_count;
        }
      }
      unsigned long long all_vertices, all_edges;
      MPI_Allreduce(&vertex_count, &all_vertices, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&edge_count, &all_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      

      // std::unordered_map<Vertex, Vertex> label_map;
      // std::vector<std::pair<Vertex, Vertex> > edge_list;
      // std::string sample_filename = result_dir + //"/" + std::to_string(0) +
      //                                      "/sampled_graph_" + std::to_string(mpi_rank);
      // std::ofstream sample(sample_filename, std::ofstream::out);
      // prunejuice::sample_graph(graph, vertex_metadata, vertex_count, label_map, edge_list, sample);
      // std::cout<<label_map.size()<<" "<<edge_list.size()<<std::endl;
      // for(auto p : label_map){
      //   sample<<"v " << p.first<<" "<<p.second<<std::endl;
      // }
      // for(auto e:edge_list){
      //   sample<<"e "<<e.first<<" "<<e.second<<std::endl;
      // }
      // return 0;

      //MPI_Barrier(MPI_COMM_WORLD); // Test
      //return 0; // Test

      /////////////////////////////////////////////////////////////////////////////

      // result
      std::string pattern_set_result_filename = result_dir + "/result_pattern_set";
      std::ofstream pattern_set_result_file;
      if (mpi_rank == 0)
      {
        pattern_set_result_file = std::ofstream(pattern_set_result_filename, std::ofstream::out);
      }
      /////////////////////////////////////////////////////////////////////////////
      // TODO: setup pattern set
      // a pattern set is a collection of directories containing pattern files

      // TODO: code indentation

      // loop over pattern set
      for (size_t ps = 0; ps < 1; ps++)
      { // TODO: for now, only reading from pattern_dir/0
        // beginning of the pattern set

        // setup pattern to search
        if (mpi_rank == 0)
        {
          std::cout << "Setting up Pattern [" << ps << "] ... " << std::endl;
        }

        std::string pattern_input_filename = pattern_dir + "/pattern";
        typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData, EdgeData> PatternGraph;
        PatternGraph pattern_graph(pattern_input, false);
        // PatternGraph pattern_graph(
        //     pattern_input_filename + "_edge",
        //     pattern_input_filename + "_vertex",
        //     pattern_input_filename + "_vertex_data",
        //     pattern_input_filename + "_edge_data",
        //     pattern_input_filename + "_stat",
        //     false, false); // TODO: improve

        // printing the pattern
        // std::cout<<"directed:"<<pattern_graph.directed<<std::endl;
        // std::cout<<"vertex_count:"<<pattern_graph.vertex_count<<std::endl;
        // std::cout<<"edge_count:"<<pattern_graph.edge_count<<std::endl;
        // std::cout<<"diameter:"<<pattern_graph.diameter<<std::endl;
        // std::cout<<"vertices:{";
        // for(auto v : pattern_graph.vertices){
        //   std::cout<<v<<", ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"vertex_degree:{";
        // for(auto v : pattern_graph.vertex_degree){
        //   std::cout<<v<<", ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"vertex_data:{";
        // for(auto v : pattern_graph.vertex_data){
        //   std::cout<<v<<", ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"edges:{";
        // for(auto v : pattern_graph.edges){
        //   std::cout<<v<<", ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"edge_ID:{";
        // for(auto v : pattern_graph.edge_ID){
        //   std::cout<<v<<", ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"edge_data:{";
        // for(auto v : pattern_graph.edge_data){
        //   std::cout<<v<<", ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"edge_list:{";
        // for(auto v : pattern_graph.edge_list){
        //   std::cout<<"("<<std::get<0>(v)<<", "<<std::get<1>(v)<<") ";
        // }
        // std::cout<<std::endl;
        // std::cout<<"neighbor count";
        // for(auto v : pattern_graph.vertex_neighbor_data_count_map){
        //   std::cout<<"{";
        //   for(auto k : v){
        //     std::cout<<"("<<k.first<<","<<k.second<<"), ";
        //   }
        //   std::cout<<"}"<<std::endl;
        // }
        // std::cout<<std::endl;

        if (mpi_rank == 0)
        {
          // std::cout<<"----- debug ----"<<std::endl;
          // for(auto v: pattern_graph.vertex_degree){
          //   std::cout<<v<<std::endl;
          // }
          std::cout << "Pattern Matching | Searching Pattern [" << ps
                    << "] : " << std::endl;
          for (auto v = 0; v < pattern_graph.vertex_count; v++)
          {
            std::cout << v << " : off-set " << pattern_graph.vertices[v]
                      << " vertex_data " << pattern_graph.vertex_data[v]
                      << " vertex_degree " << pattern_graph.vertex_degree[v] << std::endl;
            std::cout << " neighbours : ";
            for (auto e = pattern_graph.vertices[v]; e < pattern_graph.vertices[v + 1]; e++)
            {
              auto v_nbr = pattern_graph.edges[e];
              std::cout << v_nbr << ", ";
            }
            std::cout << std::endl;
            std::cout << " neighbour vertex data count : ";
            for (auto &nd : pattern_graph.vertex_neighbor_data_count_map[v])
            {
              std::cout << "(" << nd.first << ", " << nd.second << ") ";
            }
            std::cout << std::endl;
          }
          std::cout << "diameter : " << pattern_graph.diameter << std::endl;
        }
        // test print
        // std::cout<<"------------ datagraph -------------"<<std::endl;
        // for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end(); ++vitr) {
        //   vloc_type vertex = *vitr;
        //   std::cout << mpi_rank << ", l, "
        //                  << graph->locator_to_label(vertex) << ", "
        //                  << vertex_metadata[vertex] << ", "
        //                  << graph->degree(vertex) << 
        //                  ", " << graph->outgoing_degree(vertex) << ", "
        //                  << graph->incoming_degree(vertex) << std::endl;
        // }
        //MPI_Barrier(MPI_COMM_WORLD); // TODO: ?

        typedef pattern_nonlocal_constraint<Vertex, Edge, VertexData, PatternGraph>
            PatternNonlocalConstraint;

        // PatternNonlocalConstraint ptrn_util_two(pattern_graph,
        //  pattern_input_filename + "_non_local_constraints",
        //  pattern_input_filename + "vertex_non_local_constraints");

        PatternNonlocalConstraint ptrn_util_two(pattern_graph);
        // Test

        //MPI_Barrier(MPI_COMM_WORLD); // Test
        //return 0; // Test

        // initialization - per-pattern

        // initialize containers
        vertex_state_map.clear(); // important
        //vertex_rank.reset(0);
        vertex_active.reset(true); // initially all vertices are active
        //vertex_iteration.reset(0); // TODO: -1 ?
        vertex_active_edges_map.clear(); // important
        vertex_token_source_set.clear(); // clear all the sets on all the vertices
        //edge_metadata.reset(55); // Test
        //edge_active.reset(0); //  initially all edges are active / inactive

        // initialize application parameters
        bool global_init_step = true;     // TODO: Boolean
        bool global_not_finished = false; // TODO: Boolean

        bool do_nonlocal_constraint_checking = true; // TODO: Boolean

        uint64_t global_itr_count = 0;
        uint64_t active_vertices_count = 0;
        uint64_t active_edges_count = 0;
        uint64_t message_count = 0;

        std::ofstream pruned_graph_file(result_dir+"/pruned_graph_"+std::to_string(mpi_rank));

        // result
        std::string itr_result_filename = result_dir +         //"/" + std::to_string(ps) +
                                          "/result_iteration"; // TODO: improve
        std::ofstream itr_result_file(itr_result_filename, std::ofstream::out);

        std::string step_result_filename = result_dir + //"/" + std::to_string(ps)
                                           "/result_step";
        std::ofstream step_result_file(step_result_filename, std::ofstream::out);

        std::string superstep_result_filename = result_dir + //"/" + std::to_string(ps) +
                                                "/result_superstep";
        std::ofstream superstep_result_file(superstep_result_filename, std::ofstream::out);

        std::string active_vertices_count_result_filename = result_dir + //"/" + std::to_string(ps) +
                                                            "/all_ranks_active_vertices_count/active_vertices_" + std::to_string(mpi_rank);
        std::ofstream active_vertices_count_result_file(active_vertices_count_result_filename, std::ofstream::out);

        std::string active_vertices_result_filename = result_dir + //"/" + std::to_string(ps) +
                                                      "/all_ranks_active_vertices/active_vertices_" + std::to_string(mpi_rank);
        std::ofstream active_vertices_result_file(active_vertices_result_filename, std::ofstream::out);

        std::string active_edges_count_result_filename = result_dir + //"/" + std::to_string(ps) +
                                                         "/all_ranks_active_edges_count/active_edges_" + std::to_string(mpi_rank);
        std::ofstream active_edges_count_result_file(active_edges_count_result_filename, std::ofstream::out);

        std::string active_edges_result_filename = result_dir + //"/" + std::to_string(ps) +
                                                   "/all_ranks_active_edges/active_edges_" + std::to_string(mpi_rank);
        std::ofstream active_edges_result_file(active_edges_result_filename, std::ofstream::out);

        std::string paths_result_filename = result_dir + "/" +
         std::to_string(ps) + "/all_ranks_paths/paths_" + std::to_string(mpi_rank);
        std::ofstream paths_result_file(paths_result_filename, std::ofstream::out);

        std::string message_count_result_filename = result_dir +                                                //"/" + std::to_string(ps) +
                                                    "/all_ranks_messages/messages_" + std::to_string(mpi_rank); // TODO:message_count
        std::ofstream message_count_result_file(message_count_result_filename, std::ofstream::out);

        MPI_Barrier(MPI_COMM_WORLD);

        double pattern_time_start = MPI_Wtime();

        /////////////////////////////////////////////////////////////////////////////
        std::unordered_map<Vertex, Vertex> count;
        for (auto it = graph->vertices_begin(); it != graph->vertices_end(); it++){
          auto itf = count.find(vertex_metadata[*it]);
          if(itf == count.end()){
            count.insert({vertex_metadata[*it], 1});
          }else{
            itf->second ++;
          }
        }
        prunejuice::gather_map(count, 0);
        if(mpi_rank == 0){
          for(auto p : count){
            std::cout<<"("<<p.first<<","<<p.second<<"),";
          }
          std::cout<<std::endl;
        }
        // return 0;


        global_not_finished = false;

        double itr_time_start = MPI_Wtime();

        /////////////////////////////////////////////////////////////////////////////

        // mark inactive edges
        //update_edge_state();

        /////////////////////////////////////////////////////////////////////////////
        //#ifdef ENABLE_BLOCK
        // label propagation
        double label_propagation_time_start = MPI_Wtime();
        double pruning_time = 0;
#ifdef ENABLE_MY_LCC
        // --- lcc
        prunejuice::label_propagation_label_degree_filtering_bsp<
          Vertex, VertexData, VertexMetadata, BitSet, NLFMap, graph_type,
          PatternGraph, VertexStateMap, VertexActive, VertexUint8MapCollection,
          TemplateVertex>(graph, vertex_metadata, pattern_graph, vertex_state_map, vertex_active,
        vertex_active_edges_map, template_vertices, global_init_step,
        global_not_finished, superstep_result_file, message_count_result_file,
        active_vertices_count_result_file, active_edges_count_result_file);
        
        MPI_Barrier(MPI_COMM_WORLD);
        count.clear();
        for (auto & v: vertex_state_map){
          auto itf = count.find(vertex_metadata[graph->label_to_locator(v.first)]);
          if(itf == count.end()){
            count.insert({vertex_metadata[graph->label_to_locator(v.first)], 1});
          }else{
            itf->second ++;
          }
        }
        prunejuice::gather_map(count, 0);
        if(mpi_rank == 0){
          std::cout<<"after:"<<std::endl;
          for(auto p : count){
            std::cout<<"("<<p.first<<","<<p.second<<"),";
          }
          std::cout<<std::endl;
        }
#else
      prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData,
                                                           graph_type, VertexMetadata, VertexStateMap, VertexActive,
                                                           VertexUint8MapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>(graph, vertex_metadata, vertex_state_map, vertex_active,
                                                                                                                                         vertex_active_edges_map, template_vertices, pattern_graph, global_init_step,
                                                                                                                                         global_not_finished, global_itr_count, superstep_result_file,
                                                                                                                                         active_vertices_count_result_file, active_edges_count_result_file, message_count_result_file);    
#endif

        MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
        double label_propagation_time_end = MPI_Wtime();
        unsigned long long LCC_remain_vertices, LCC_remain_edges;

        get_pruned_graph_size(graph, vertex_state_map, vertex_active_edges_map, LCC_remain_vertices, LCC_remain_edges);

        if(mpi_rank==0){
          std::cout<<"LCC time and size:"<<label_propagation_time_end-label_propagation_time_start<<":"<<LCC_remain_vertices<<":"<<LCC_remain_edges<<std::endl;
          pruning_time += label_propagation_time_end-label_propagation_time_start;
        }


        //#endif

        /////////////////////////////////////////////////////////////////////////////

        if (global_init_step)
        { // Important
          global_init_step = false;
        }

        global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

        // global verification - are all vertex_state_maps empty
        // false - no active vertex left, true - active vertices left
        bool global_active_vertex = true; //vertex_state_map.size() < 1 ? false : true;

        // forced token passing
        if (global_itr_count == 0)
        {
          global_not_finished = true; // TODO: for load balancing experiments, ?
        }

        // Test
        // toekn passing
        double token_passing_time_start = MPI_Wtime();

        if (ptrn_util_two.input_patterns.size() < 1)
        {
          do_nonlocal_constraint_checking = false;
        }
        // std::cout<<"-----> non constraint size:"<<ptrn_util_two.input_patterns.size()<<" "<<global_not_finished<<std::endl;
        if (do_nonlocal_constraint_checking && global_not_finished)
        { // TODO: do we need this?

          global_not_finished = false;

          VertexUint8Map token_source_map; // per vertex state

          VectorBoolean pattern_found(ptrn_util_two.input_patterns.size(), 0);              // per rank state
          VectorBoolean pattern_token_source_found(ptrn_util_two.input_patterns.size(), 0); // per rank state

          // TODO: code indentation
          // loop over the constraints and run token passing
          for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++)
          {

            // TODO: only output subgraphs when doing enumeration
            // result
            std::string paths_result_filename = result_dir + //"/" +
                                                //std::to_string(ps) + "/all_ranks_paths/paths_" +
                                                //std::to_string(ps) +
                                                "/all_ranks_subgraphs/subgraphs_" +
                                                std::to_string(pl) + "_" + std::to_string(mpi_rank);
            std::ofstream paths_result_file(paths_result_filename, std::ofstream::out);

            bool token_source_deleted = false;

            // TODO: This is actually a bad design. It should be one object per entry.
            // ptrn_util_twois the object

            // setup pattern
            auto pattern_tp = std::get<0>(ptrn_util_two.input_patterns[pl]);
            auto pattern_indices_tp = std::get<1>(ptrn_util_two.input_patterns[pl]);
            auto pattern_cycle_length_tp = std::get<2>(ptrn_util_two.input_patterns[pl]); // uint
            auto pattern_valid_cycle_tp = std::get<3>(ptrn_util_two.input_patterns[pl]);  // boolean
                                                                                          //  auto pattern_interleave_label_propagation_tp = std::get<4>(ptrn_util_two.input_patterns[pl]); // boolean
                                                                                          //--  auto pattern_seleted_edges_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean
                                                                                          //  auto pattern_selected_vertices_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean

            auto pattern_selected_vertices_tp = 0; // TODO: remove

            auto pattern_is_tds_tp = std::get<4>(ptrn_util_two.input_patterns[pl]);                       // boolean
            auto pattern_interleave_label_propagation_tp = std::get<5>(ptrn_util_two.input_patterns[pl]); // boolean

            auto pattern_enumeration_tp = ptrn_util_two.enumeration_patterns[pl];
            auto pattern_aggregation_steps_tp = ptrn_util_two.aggregation_steps[pl];

            // TODO: read from file / remove
            auto pattern_selected_edges_tp = false;     // boolean
            auto pattern_mark_join_vertex_tp = false;   // boolean
            auto pattern_ignore_join_vertex_tp = false; // boolean
            size_t pattern_join_vertex_tp = 0;          // TODO:
            //bool do_tds_tp = false;

            message_count = 0;

            if (mpi_rank == 0)
            {
              std::cout << "Token Passing [" << pl << "] | Searching Subpattern : ";
              //pattern_util<VertexData>::output_pattern(pattern_tp);
              PatternNonlocalConstraint::output_pattern(pattern_tp);
              std::cout << "Token Passing [" << pl << "] | Vertices : ";
              //pattern_util<VertexData>::output_pattern(pattern_indices_tp);
              PatternNonlocalConstraint::output_pattern(pattern_indices_tp);
              std::cout << "Token Passing [" << pl << "] | Arguments : "
                        << pattern_cycle_length_tp << " "
                        << pattern_valid_cycle_tp << " "
                        << pattern_interleave_label_propagation_tp << " "
                        //--      << pattern_seleted_edges_tp << std::endl; // Test
                        << pattern_selected_vertices_tp << std::endl; // Test

              std::cout << "Token Passing [" << pl << "] | Enumeration Indices : ";
              //pattern_util<Vertex>::output_pattern(pattern_enumeration_tp);
              PatternNonlocalConstraint::output_pattern(pattern_enumeration_tp);
              std::cout << "Token Passing [" << pl << "] | Agreegation Steps : TODO" << std::endl;
              for(auto p : pattern_aggregation_steps_tp){
                std::cout<<(int)p<<":";
              }
              std::cout<<std::endl;
              // std::cout<<"AGGREGATION:"<<pattern_aggregation_steps_tp[0]<<std::endl;
              //PatternNonlocalConstraint::output_pattern(pattern_aggregation_steps_tp); // TODO:
            }

            //return 0; // Test

            // initialize containers
            unsigned long long local_active=0, global_active;
            if (!pattern_selected_vertices_tp)
            {
              token_source_map.clear();        // Important
              vertex_token_source_set.clear(); // Important
            }
            else
            {

              token_source_map.clear(); // Important

              for (vitr_type vitr = graph->vertices_begin();
                   vitr != graph->vertices_end(); ++vitr)
              {
                local_active ++;
                vloc_type vertex = *vitr;
                if (vertex_active[vertex] &&
                    (vertex_metadata[vertex] == pattern_tp[pattern_tp.size() - 1]))
                {
                  //std::cout << mpi_rank << " " << graph->locator_to_label(vertex) << " [l] " << vertex_token_source_set[vertex].size() << std::endl; // Test
                  continue;
                }
                else
                {
                  vertex_token_source_set[vertex].clear();
                }
              }

              for (vitr_type vitr = graph->delegate_vertices_begin();
                   vitr != graph->delegate_vertices_end(); ++vitr)
              {
                vloc_type vertex = *vitr;
                if (vertex_active[vertex] &&
                    (vertex_metadata[vertex] == pattern_tp[pattern_tp.size() - 1]))
                {
                  //std::cout << mpi_rank << " " << graph->locator_to_label(vertex) << " [d] " << vertex_token_source_set[vertex].size() << std::endl; // Test
                  continue;
                }
                else
                {
                  vertex_token_source_set[vertex].clear();
                }
              }

            } // if pattern_selected_vertices_tp
            

            local_active = 0;
            for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end(); ++vitr){
              vloc_type vertex = *vitr;
              if (vertex_active[vertex])
              {
                local_active ++;
                //std::cout << mpi_rank << " " << graph->locator_to_label(vertex) << " [l] " << vertex_token_source_set[vertex].size() << std::endl; // Test
                continue;
              }
            }
            MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            //--  //token_source_edge_set.clear(); // Important // edge aware

            time_start = MPI_Wtime();

#ifdef TP_ASYNC
            if (pattern_is_tds_tp)
            {
#ifdef ENABLE_TDS
              MPI_Barrier(MPI_COMM_WORLD);
              double start_tds_time = MPI_Wtime();
#endif
              prunejuice::token_passing_pattern_matching<graph_type, Vertex, Edge, VertexData,
                                                         EdgeData, VertexMetadata, EdgeMetadata, VertexActive,
                                                         VertexUint8MapCollection, TemplateVertex, VertexStateMap, PatternGraph,
                                                         /*PatternUtilities*/ PatternNonlocalConstraint, VertexUint8Map, VertexSetCollection,
                                                         DelegateGraphVertexDataSTDAllocator, Boolean, BitSet>(graph, vertex_metadata, vertex_active, vertex_active_edges_map,
                                                                                                               template_vertices, vertex_state_map, pattern_graph, ptrn_util_two, pl,
                                                                                                               token_source_map, vertex_token_source_set,
                                                                                                               pattern_found[pl], tp_vertex_batch_size, paths_result_file, message_count);

#ifdef ENABLE_TDS
              MPI_Barrier(MPI_COMM_WORLD);
              if(mpi_rank == 0){
                  std::cout<<"TDS search time:"<<MPI_Wtime()-start_tds_time<<std::endl;
              }
              return 0;
#endif

            }
            else
            {
              // prunejuice::<graph_type, VertexMetadata, decltype(pattern_tp), decltype(pattern_indices_tp), uint8_t, PatternGraph,
              //                                            VertexStateMap, VertexUint8Map, edge_data_t,
              //                                            VertexSetCollection, VertexActive, TemplateVertex, VertexUint8MapCollection, BitSet>(graph, vertex_metadata, pattern_tp,
              //                                                                                                                                 pattern_indices_tp, vertex_rank, pattern_graph, vertex_state_map,
              //                                                                                                                                 token_source_map, pattern_cycle_length_tp, pattern_valid_cycle_tp,
              //                                                                                                                                 pattern_found[pl], *edge_data_ptr, vertex_token_source_set, vertex_active,
              //                                                                                                                                 template_vertices, vertex_active_edges_map, pattern_selected_vertices_tp, //);
              //                                                                                                                                 pattern_selected_edges_tp, pattern_mark_join_vertex_tp,
              //                                                                                                                                 pattern_ignore_join_vertex_tp, pattern_join_vertex_tp, message_count);

              prunejuice::token_passing_pattern_filtering<graph_type, Vertex, Edge, VertexData,
                                                         EdgeData, VertexMetadata, EdgeMetadata, VertexActive,
                                                         VertexUint8MapCollection, TemplateVertex, VertexStateMap, PatternGraph,
                                                         /*PatternUtilities*/ PatternNonlocalConstraint, VertexUint8Map, VertexSetCollection,
                                                         DelegateGraphVertexDataSTDAllocator, Boolean, BitSet>(graph, vertex_metadata, vertex_active, vertex_active_edges_map,
                                                                                                               template_vertices, vertex_state_map, pattern_graph, ptrn_util_two, pl,
                                                                                                               token_source_map, vertex_token_source_set,
                                                                                                               pattern_found[pl], tp_vertex_batch_size, paths_result_file, message_count);

            }
#endif
            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
            time_end = MPI_Wtime();

            unsigned long long NLCC_remain_vertices, NLCC_remain_edges;
            get_pruned_graph_size(graph, vertex_state_map, vertex_active_edges_map, NLCC_remain_vertices, NLCC_remain_edges);

            if (mpi_rank == 0)
            {
              std::cout << "Pattern Matching Time | Token Passing (Traversal) ["
                        << pl << "] : " << time_end - time_start << std::endl;
              std::cout<<"NLCC-"<<pl<<":"<<time_end-time_start<<":"<<NLCC_remain_vertices<<":"<<NLCC_remain_edges<<std::endl;
              pruning_time += time_end-time_start;
            }



            uint64_t remove_count = 0;
            size_t token_source_pattern_indices_tp_index = 0; // TODO: this is confusing, update

            bool is_token_source_map_not_empty = havoqgt::mpi_all_reduce(!token_source_map.empty(), std::greater<uint8_t>(), MPI_COMM_WORLD); // TODO: less did not work?
            MPI_Barrier(MPI_COMM_WORLD);  
#ifdef TP_ASYNC

            // remove the invalid (token source) vertices from the vertex_state_map
            // for delegates, set vertex_active to false

            // TODO: In the case, a vertex is on multiple cycles/chains (not as the token source)
            // only invalidate it as a token source, but do not remove it from the vertex_state_map

                                                                                                                // TODO: might not need this here
            unsigned long long token_size = 0;
            unsigned long long local_token_size = token_source_map.size();
            MPI_Allreduce(&local_token_size, &token_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            
            pattern_token_source_found[pl] = is_token_source_map_not_empty;
            
            // std::string token_file_name = result_output+"/token_"+std::to_string(mpi_rank);
            // std::ofstream token_file(token_file_name, std::ofstream::out);
            local_active = 0;
            for (auto &s : token_source_map){
              
              local_active ++;
              if (!s.second){
                auto v_locator = graph->label_to_locator(s.first);
                BitSet v_template_vertices(template_vertices[v_locator]);
                if (v_template_vertices.none())
                {
                  continue;
                }
                // token_file<<mpi_rank<<":"<<s.first<<":"<<(int)s.second<<":"<<v_template_vertices<<":";
                //pattern_indices_tp[0]; // token source template vertex ID

                if (v_template_vertices.test(pattern_indices_tp[token_source_pattern_indices_tp_index]))
                {
                  assert(pattern_indices_tp[token_source_pattern_indices_tp_index] < max_template_vertex_count); // Test
                  v_template_vertices.reset(pattern_indices_tp[token_source_pattern_indices_tp_index]);
                  template_vertices[v_locator] = v_template_vertices.to_ulong();
                }
                // token_file<<v_template_vertices<<":"<<v_template_vertices.none()<<std::endl;
                if (v_template_vertices.none())
                {
                  // token_file<<"-------"<<std::endl;
                  vertex_active[graph->label_to_locator(s.first)] = false;
                }

                if (!global_not_finished)
                {
                  global_not_finished = true;
                }

                if (!token_source_deleted)
                {
                  token_source_deleted = true;
                }
              }
            }
            
            get_pruned_graph_size(graph, vertex_state_map, vertex_active_edges_map, NLCC_remain_vertices, NLCC_remain_edges);
            
            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

            vertex_active.all_min_reduce();
            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?

            // TODO: this is a temporary patch, forcing all the delegates to have no identity
            for (vitr_type vitr = graph->delegate_vertices_begin();
                 vitr != graph->delegate_vertices_end(); ++vitr)
            {
              auto vertex = *vitr;
              if (vertex.is_delegate() && (graph->master(vertex) == mpi_rank))
              {
                continue; // skip the controller
              }
              else
              {
                auto find_vertex = vertex_state_map.find(graph->locator_to_label(vertex));
                if (find_vertex == vertex_state_map.end())
                {
                  template_vertices[vertex] = 0;
                }
              }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            get_pruned_graph_size(graph, vertex_state_map, vertex_active_edges_map, NLCC_remain_vertices, NLCC_remain_edges);
            
            template_vertices.all_max_reduce(); // ensure all the delegates have the same value as the controller
            MPI_Barrier(MPI_COMM_WORLD);        // TODO: do we need this here?
            // std::string result_active_name = result_output+"/active_"+std::to_string(mpi_rank);
            // std::ofstream active(result_active_name, std::ofstream::out); 

            // for(auto vitr=graph->vertices_begin(); vitr!= graph->vertices_end(); ++ vitr){
            //     if(vertex_active[*vitr]){
            //         active<<graph->locator_to_label(*vitr)<<std::endl;
            //     }
            // }
            // MPI_Allreduce(&local_active, &global_active, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            // if(mpi_rank == 0){
            //     std::cout<<"active size:"<<global_active<<std::endl;
            // }
            // remove from vertex_state_map
            for (auto &s : token_source_map)
            {
              auto v_locator = graph->label_to_locator(s.first);
              if (!vertex_active[v_locator])
              {
                auto find_vertex = vertex_state_map.find(s.first);

                if (find_vertex != vertex_state_map.end())
                {
                  vertex_state_map.erase(find_vertex);
                  remove_count++;
                  // if (vertex_state_map.erase(s.first) < 1)
                  // { // s.first is the vertex
                  //   std::cerr << "Error: failed to remove an element from the map."
                  //             << std::endl;
                  // }
                  // else
                  // {
                  //   remove_count++;
                  //   //if (!global_not_finished) {
                  //   //  global_not_finished = true;
                  //   //}
                  // }
                }
              }
            }

#endif // ifdef TP_ASYNC
            // return 0;
            // Important : This may slow down things -only for presenting results
            active_vertices_count = 0;
            active_edges_count = 0;

            for (auto &v : vertex_state_map)
            {
              auto v_locator = graph->label_to_locator(v.first);
              if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank))
              {
                active_vertices_count++;

                // edges
                active_edges_count += vertex_active_edges_map[v_locator].size();
              }
              else if (!v_locator.is_delegate())
              {
                active_vertices_count++;

                // edges
                active_edges_count += vertex_active_edges_map[v_locator].size();
              }
            }

            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

            if (is_token_source_map_not_empty)
            {
              havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
              if (mpi_rank == 0)
              {
                //    for (size_t pl = 0; pl < ptrn_util_two.input_patterns.size(); pl++) {
                std::string s = pattern_found[pl] == 1 ? "True" : "False";
                std::cout << "Token Passing [" << pl << "] | Found Subpattern : " << s << std::endl;
                //    }
              }

              MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

              // verify global token source deleted status
              //token_source_deleted = havoqgt::mpi::mpi_all_reduce(token_source_deleted, std::logical_or<bool>(), MPI_COMM_WORLD); // does not work
              token_source_deleted = havoqgt::mpi_all_reduce(token_source_deleted, std::greater<uint8_t>(), MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
              if (mpi_rank == 0)
              {
                std::cout << "Token Passing [" << pl << "] | Token Source Deleted Status : ";
                if (token_source_deleted)
                {
                  std::cout << "Deleted" << std::endl;
                }
                else
                {
                  std::cout << "Not Deleted" << std::endl;
                }
              }
            }
            else
            {
              if (mpi_rank == 0)
              {
                std::cout << "Token Passing [" << pl << "] | No Token Source Found" << std::endl;
              }
            } // is_token_source_map_not_empty

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

            // interleave token passing with label propagation
            if (token_source_deleted && pattern_interleave_label_propagation_tp)
            {

              bool global_not_finished_dummy = true; // TODO: do we need this?

              // lable propagation
              label_propagation_time_start = MPI_Wtime();

// debug
            get_pruned_graph_size(graph, vertex_state_map, vertex_active_edges_map, NLCC_remain_vertices, NLCC_remain_edges);

#ifdef LCC_OLD
              prunejuice::label_propagation_pattern_matching_bsp<graph_type, VertexMetadata, VertexData, decltype(pattern), decltype(pattern_indices),
                                                                 /*VertexRank*/ uint8_t, VertexActive, /*VertexIteration*/ uint8_t, VertexStateMap, PatternGraph, BitSet, TemplateVertex,
                                                                 VertexUint8MapCollection>(graph, vertex_metadata, pattern, pattern_indices, vertex_rank, vertex_active,
                                                                                           vertex_iteration, vertex_state_map, pattern_graph, global_init_step, global_not_finished,
                                                                                           global_itr_count, superstep_result_file, active_vertices_count_result_file, active_edges_count_result_file,
                                                                                           template_vertices, vertex_active_edges_map, message_count_result_file);
#endif

#ifdef ENABLE_MY_LCC
        prunejuice::label_propagation_label_degree_filtering_bsp<
          Vertex, VertexData, VertexMetadata, BitSet, NLFMap, graph_type,
          PatternGraph, VertexStateMap, VertexActive, VertexUint8MapCollection,
          TemplateVertex>(graph, vertex_metadata, pattern_graph, vertex_state_map, vertex_active,
        vertex_active_edges_map, template_vertices, global_init_step,
        global_not_finished, superstep_result_file, message_count_result_file,
        active_vertices_count_result_file, active_edges_count_result_file);
#else
              prunejuice::label_propagation_pattern_matching_bsp<Vertex, VertexData,
                                                                 graph_type, VertexMetadata, VertexStateMap, VertexActive,
                                                                 VertexUint8MapCollection, TemplateVertexBitSet, TemplateVertex, PatternGraph>(graph, vertex_metadata, vertex_state_map, vertex_active,
                                                                                                                                               vertex_active_edges_map, template_vertices, pattern_graph, global_init_step,
                                                                                                                                               global_not_finished, global_itr_count, superstep_result_file,
                                                                                                                                               active_vertices_count_result_file, active_edges_count_result_file,
                                                                                                                                               message_count_result_file);

#endif
              MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here
              label_propagation_time_end = MPI_Wtime();
              get_pruned_graph_size(graph, vertex_state_map, vertex_active_edges_map, NLCC_remain_vertices, NLCC_remain_edges);
            }

            // result
            paths_result_file.close();

            MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here? // New

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//

          } // for - loop over the constraints and run token passing

          // pattren found ? // TODO: write output to file
          ///havoqgt::mpi::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
          havoqgt::mpi_all_reduce_inplace(pattern_found, std::greater<uint8_t>(), MPI_COMM_WORLD);
          MPI_Barrier(MPI_COMM_WORLD); // TODO: do we need this here?
        }
        else
        {
          global_not_finished = false;
        } // do token passing ?

        global_not_finished = havoqgt::mpi_all_reduce(global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD); // TODO: might not need this here

        // if (mpi_rank == 0)
        // {
        //   std::cout << "Pattern Matching | Global Finished Status : ";
        //   if (global_not_finished)
        //   {
        //     std::cout << "Continue" << std::endl;
        //   }
        //   else
        //   {
        //     std::cout << "Stop" << std::endl;
        //   }
        // }
        // finished
        // std::cout<<"----------------------------------------"<<vertex_state_map.size()<<std::endl;
        // active_vertices_count = 0;
        // for (auto& v : vertex_state_map) {
        //   auto v_locator = graph->label_to_locator(v.first);
        //   if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
        //     active_vertices_count++;
        //     std::cout<<"vertex:"<<v.first<<std::endl;
        //   } else if (!v_locator.is_delegate()) {
        //     active_vertices_count++;
        //     std::cout<<"vertex:"<<v.first<<std::endl;
        //   }
        // }
      //   std::cout << "------------ pattern graph -------------" << std::endl;
      // std::cout << "Directed: " << pattern_graph.directed << std::endl;
      // for (auto v = 0; v < pattern_graph.vertex_count; v++) {
      //   std::cout << v << ": off-set " << pattern_graph.vertices[v]
      //             << " vertex_data " << pattern_graph.vertex_data[v]
      //             << " vertex_degree " << pattern_graph.vertex_degree[v]
      //             // << " vertex_outgoing_degree " << pattern_graph.vertex_degree[v]
      //             // << " vertex_incoming_degree " << pattern_graph.vertex_incoming_degree[v]
      //             << std::endl;
      //   std::cout << " neighbors: ";
      //   for (auto e = pattern_graph.vertices[v];
      //        e < pattern_graph.vertices[v + 1]; e++) {
      //     auto v_nbr = pattern_graph.edges[e];
      //     std::cout << v_nbr << ", ";
      //   }
      //   std::cout << std::endl;
      // }

        /////////////////////////////////////////////////////////////////////////////

        double itr_time_end = MPI_Wtime();
        //  if(mpi_rank == 0) { //TODO: sum of LCC and NLCC iterations
        //    std::cout << "Pattern Matching Time | Pattern [" << ps
        //      << "] | Iteration [" << global_itr_count << "] : "
        //      << itr_time_end - itr_time_start << std::endl;
        //  }

        // result
        if (mpi_rank == 0)
        {
          // iteration number, time
          itr_result_file << global_itr_count << ", "
                          << (itr_time_end - itr_time_start) << "\n";
        }

        global_itr_count++; //TODO: sum of LCC and NLCC iterations

        // global termination
        // Test
        //if (global_itr_count > 0) {
        //  global_not_finished = false;
        //}
        //MPI_Barrier(MPI_COMM_WORLD);
        // Test

        ///  } while (global_not_finished); // application loop

        /////////////////////////////////////////////////////////////////////////////

        MPI_Barrier(MPI_COMM_WORLD);
        double pattern_time_end = MPI_Wtime();
        // if (mpi_rank == 0)
        // {
        //   std::cout << "Pattern Matching Time | Pattern [" << ps << "] : "
        //             << pattern_time_end - pattern_time_start << std::endl;
        // }
        
        // std::cout<<"------------- local prune ---------"<<std::endl;
        // for (vitr_type vitr = graph->vertices_begin(); vitr != graph->vertices_end(); ++vitr) {
        //   vloc_type vertex = *vitr;
        //   if (vertex_active[vertex]) {
        //     std::cout << mpi_rank << ", l, " << graph->locator_to_label(vertex) << ", " << graph->degree(vertex) << ", " << vertex_metadata[vertex] << "\n";  
        // //  vertex_active_count++;
        //   } else { 
        //   // vertex_inactive_count++;
        //   }  
        // }


        /*
        * starting enumeration process
        */
        vector<Vertex> ordering;
        generate_ordering_vf3(pattern_graph, graph, vertex_metadata, mpi_rank, mpi_size, ordering);
        if(mpi_rank == 0){
          std::cout<<"------------->ordering:"<<std::endl;
          for(auto v : ordering){
            std::cout<<v<<" ";
          }
          std::cout<<std::endl;
        }

        std::unordered_map<Vertex, std::set<Vertex>> active_map;
        for(auto &v : vertex_state_map){
          auto v_locator = graph->label_to_locator(v.first);
          active_map.insert({v.first, {}});
          for(auto &n : vertex_active_edges_map[v_locator]){
            active_map[v.first].insert(n.first);
          }
        }
        // for(auto loc : vertex_active_edges_map){
        //   if(loc.owner() == mpi_rank){
        //     Vertex s = graph.locator_to_label(loc);
        //     active_map.insert(s, {});
        //     for(auto d : vertex_active_edges_map[loc]){
        //       active_map[s].insert(d);
        //     }
        //   }
        // }
        
        std::string result_file_name = result_output+"/pattern_result_"+std::to_string(mpi_rank);
        std::string pruned_file_name = result_output+"/pruned_graph_"+std::to_string(mpi_rank);
        std::ofstream result(result_file_name, std::ofstream::out);
        std::ofstream pruned_graph(pruned_file_name, std::ofstream::out);
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef ENABLE_CACHE
        Query_cache qcache(vertex_state_map, vertex_active_edges_map, graph, mpi_rank, mpi_size, all_vertices);
#endif
        double start_enumrtate_time = MPI_Wtime();
        if(mpi_rank == 0){
          std::cout<<"------- start enumeration -------"<<std::endl;
        }
        unsigned long long result_count = 0;
        std::vector<unsigned long long> local_computation{0, 0, 0};
        std::vector<double> time_col{0,0,0,0,0};
#ifdef ENABLE_CACHE
        pattern_enumeration(graph, &pattern_graph, ordering, vertex_metadata, vertex_state_map, template_vertices, vertex_state_map, vertex_active_edges_map, result, result_count, active_map, time_col, local_computation, qcache, pruned_graph); //, local_computation
#else
        pattern_enumeration(graph, &pattern_graph, ordering, vertex_metadata, vertex_state_map, template_vertices, vertex_state_map, vertex_active_edges_map, result, result_count, active_map, time_col, local_computation, pruned_graph); //, local_computation
#endif
        MPI_Barrier(MPI_COMM_WORLD);
        double end_enumrtate_time = MPI_Wtime();
        unsigned long long all_result_count;
        MPI_Allreduce(&result_count, &all_result_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if(mpi_rank == 0){
          std::cout<<"enumeration_time:"<<end_enumrtate_time-start_enumrtate_time<<std::endl;
          std::cout<<"total result counts:"<<all_result_count<<std::endl;
          std::cout<<"total_pruning_time:"<<pruning_time<<std::endl;
        }
        // std::vector<unsigned long long> recbuf;
        // recbuf.resize(mpi_size*3);
        // MPI_Gather(&local_computation, 3, MPI_UNSIGNED_LONG_LONG, &recbuf[0], 3, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        // for(int i=0;i<time_col.size();++i){
        //   double times;
        //   // std::cout<<"time_col["<<i<<"]="<<time_col[i]<<std::endl;
        //   MPI_Allreduce(&time_col[i], &times, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //   if(mpi_rank == 0){
        //     std::cout<<"time["<<i<<"]:"<<times<<std::endl;
        //   }
        // }

        /////////////////////////////////////////////////////////////////////////////
        unsigned long long remain_vertices=0,remain_edges=0;
#ifdef PRINT_PRUNED_GRAPH
        for (auto &v : vertex_state_map)
        {
          auto v_locator = graph->label_to_locator(v.first);
          BitSet v_template_vertices(template_vertices[v_locator]);
          pruned_graph_file<<"v "<<v.first<<" "<<vertex_metadata[v_locator]<<" "<<v_template_vertices<<std::endl;
          ++ remain_vertices;
          // edges
          for (auto &n : vertex_active_edges_map[v_locator])
          {
            pruned_graph_file<<"e "<<v.first<<" "<<n.first<<std::endl;
            ++ remain_edges;
          }
        }
#endif
        unsigned long long remain_all_vertices, remain_all_edges;
        MPI_Allreduce(&remain_vertices, &remain_all_vertices, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&remain_edges, &remain_all_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if(mpi_rank == 0){
          std::cout<<"all_vertices:"<<all_vertices<<"|all_edges:"<<all_edges/2<<std::endl;
          std::cout<<"remain_vertices:"<<remain_all_vertices<<"|remain_edges:"<<remain_all_edges/2<<std::endl;
        }

        // metrics about query cache
        unsigned long long all_states, redundant_states, extra_states;
        MPI_Allreduce(&(local_computation[0]), &all_states, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(local_computation[1]), &redundant_states, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(local_computation[2]), &extra_states, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if(mpi_rank == 0){
          std::cout<<"all_visitors:"<<all_states<<std::endl;
          std::cout<<"redundant_states:"<<redundant_states<<std::endl;
          std::cout<<"extra_states:"<<extra_states<<std::endl;
        }
#ifdef ENABLE_CACHE
        unsigned long long total_q, hit_q;
        // std::cout<<"-->"<<mpi_rank<<":"<<qcache.total_count<<":"<<qcache.hit_count<<std::endl;
        MPI_Allreduce(&(qcache.total_count), &total_q, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&(qcache.hit_count), &hit_q, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if(mpi_rank == 0){
          std::cout<<"hit_rate:"<<(float)hit_q/total_q<<":"<<hit_q<<":"<<total_q<<std::endl;
        }
        
#endif

        pruned_graph_file.close();

        // end of an elemnet in the pattern set
      } // for - loop over pattern set

      // if (mpi_rank == 0)
      // {
      //   pattern_set_result_file.close(); // close file
      // }

    } // pattern matching

  } // havoqgt_init
  ;
  // END Main MPI
  ///havoqgt::havoqgt_finalize();

  return 0;
} // main
