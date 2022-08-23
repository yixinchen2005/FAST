#include <havoqgt/cache_utilities.hpp>
#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
//#include <havoqgt/gen_preferential_attachment_edge_list.hpp>

#include <prunejuice/vf3.hpp>

#include <prunejuice/algorithm_state.hpp>
//#include <prunejuice/local_constraint_checking.hpp>
//#include <prunejuice/non_local_constraint.hpp>
//#include <prunejuice/non_local_constraint_checking_tds_batch.hpp>
//#include <prunejuice/non_local_constraint_checking_unique.hpp>
#include <prunejuice/template.hpp>

//#include <assert.h>

//#include <algorithm>
#include <bitset>
//#include <deque>
//#include <fstream>
//#include <functional>
//#include <iostream>
//#include <memory>
//#include <set>
//#include <sstream>
#include <string>
//#include <unordered_map>
//#include <utility>
//#include <vector>

//#include <boost/bimap.hpp>
//#include <boost/bind.hpp>
//#include <boost/function.hpp>
//#include <boost/interprocess/managed_heap_memory.hpp>

#include <metadata/vertex_data_db.hpp>
#include <metadata/vertex_data_db_degree.hpp>

using namespace havoqgt;

typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

/******************** functions for parsing the command line
 * *************************/
void usage() {
  if (comm_world().rank() == 0) {
    std::cerr
        << "Usage: -i <string> -p <string> -o <string>\n"
        << " -i <string>   - input graph base filename (required)\n"
        << " -b <string>   - backup graph base filename. If set, \"input\" "
           "graph will be deleted if it exists\n"
        << " -v <string>   - vertex metadata base filename (optional, Default "
           "is degree based metadata)\n"
        << " -e <string>   - edge metadata base filename (optional)\n"
        << " -p <string>   - pattern base directory (required)\n"
        << " -o <string>   - output base directory (required)\n"
        << " -x <int>      - Token Passing batch count (optional, Default/min "
           "batch count is 1, max batch count is "
        << comm_world().size() << "\n"
        << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_input,
                    std::string& backup_graph_input,
                    std::string& vertex_metadata_input,
                    std::string& edge_metadata_input,
                    std::string& pattern_input, std::string& result_output,
                    uint64_t& tp_vertex_batch_size) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD Line :";
    for (int i = 0; i < argc; ++i) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  bool           print_help = false;
  std::bitset<3> required_input;
  required_input.reset();

  int c;

  while ((c = getopt(argc, argv, "i:b:v:e:p:o:x:h ")) != EOF) {
    switch (c) {
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
        if (tp_vertex_batch_size < 1 ||
            tp_vertex_batch_size > comm_world().size()) {
          print_help = true;
        } else if (tp_vertex_batch_size > 1) {
          tp_vertex_batch_size = comm_world().size() / tp_vertex_batch_size;
        } else {
          tp_vertex_batch_size = comm_world().size();
        }
        break;
      case '?':
        std::cerr << "Unrecognized Option : " << c << ", Ignore." << std::endl;
        print_help = true;
        break;
      default:
        abort();
    }
  }

  if (print_help || !required_input.all()) {
    usage();
    exit(-1);
  }
}

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);

  // this already contains the basic initializations for a MPI program
  havoqgt::init(&argc, &argv);

  {
    // initialize $mpi_rank and $mpi_size for each MPI process
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    // output the information of the MPI program
    if (mpi_rank == 0) {
      std::cout << "MPI Initialized With " << mpi_size << " Ranks."
                << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // parse the command line
    std::string graph_input;
    std::string backup_graph_input;
    std::string vertex_metadata_input;
    std::string edge_metadata_input;
    std::string pattern_input;
    std::string result_output;
    uint64_t    tp_vertex_batch_size = comm_world().size();
    parse_cmd_line(argc, argv, graph_input, backup_graph_input,
                   vertex_metadata_input, edge_metadata_input, pattern_input,
                   result_output, tp_vertex_batch_size);

    std::string pattern_dir = pattern_input;
    std::string result_dir  = result_output;

    // debug
    std::string debug_file_name =
        result_dir + "/debug_" + std::to_string(mpi_rank);
    std::ofstream debug_shell(debug_file_name, std::ofstream::out);

    MPI_Barrier(MPI_COMM_WORLD);

    /******************** load the graph ***********************************/
    if (mpi_rank == 0) {
      std::cout << "Loading Graph ... " << std::endl;
    }
    // the name of the input graph must not be empty
    if (backup_graph_input.size() > 0) {
      // this is used for creating a distributed graph
      distributed_db::transfer(backup_graph_input.c_str(), graph_input.c_str());
    }
    // and this is to open the input graph $ddb
    distributed_db ddb(db_open_read_only(), graph_input.c_str());

    // using pointer $graph_type to access the data graph $ddb
    auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

    // TODO: figure out a way to get it from graph_type
    // see edge_data_value_type in parallel_edge_list_reader.hpp
    typedef uint8_t edge_data_type;
    typedef graph_type::edge_data<edge_data_type,
                                  distributed_db::allocator<edge_data_type>>
         edge_data_t;
    auto edge_data_ptr =
        ddb.get_manager()->find<edge_data_t>("graph_edge_data_obj").first;

    MPI_Barrier(MPI_COMM_WORLD);
    if (mpi_rank == 0) {
      std::cout << "Done Loading Graph." << std::endl;
    }

    /******************* isomorphism matching ************************/
    /******* basic typedef *************/
    typedef uint64_t Vertex;
    typedef uint64_t Edge;
    typedef uint64_t VertexData;
    typedef uint64_t EdgeData;
    typedef uint64_t VertexClass;
    typedef graph_type::vertex_data<VertexData, std::allocator<VertexData>>
        VertexMetadata;

    VertexMetadata vertex_metadata(
        *graph);  // the data in the array are separated in different
                  // process,thus
    // each process can only access the metadata they belong to
    MPI_Barrier(MPI_COMM_WORLD);

    /******************* construct the vertex data db ****************/
    if (vertex_metadata_input.size() > 0) {
      vertex_data_db_nostdfs<graph_type, VertexMetadata, Vertex, VertexData>(
          graph, vertex_metadata, vertex_metadata_input, 10000);
    } else {
      vertex_data_db_degree<graph_type, VertexMetadata, Vertex, VertexData>(
          graph, vertex_metadata);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /********* precalculation ************/
    // read pattern graph
    std::string pattern_input_basename = pattern_dir + "/pattern";
    typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData, EdgeData>
                 PatternGraph;
    PatternGraph pattern_graph(
        pattern_input_basename + "_edge", pattern_input_basename + "_vertex",
        pattern_input_basename + "_vertex_data",
        pattern_input_basename + "_edge_data", pattern_input_basename + "_stat",
        false, false);  // TODO: improve
    // precalculation can be done in a serial way for its low time comsuption
    std::unordered_map<Vertex, double>          Pf;
    std::unordered_map<VertexData, VertexClass> classMap;  //[label -> class]
    std::unordered_map<Vertex, Vertex>
        parentMap;  //[vertex -> parent of the vertex]
    std::unordered_map<int, std::unordered_map<VertexClass, std::set<Vertex>>>
        P1set;
    // because of the graph being undirected, $S1set is depricated
    std::unordered_map<int, std::unordered_map<VertexClass, std::set<Vertex>>>
                   S1set;
    std::vector<Vertex> NodeSequence;

    MPI_Barrier(MPI_COMM_WORLD);
    double time_start = MPI_Wtime();
    computeProbabilities(pattern_graph, graph, Pf, mpi_rank, mpi_size,
                              vertex_metadata, classMap);
    generateNodeSequence(Pf, pattern_graph, NodeSequence);
    /* setting the classMap
     * TO DO:
     * by default each label maps to a unique class
     * examples: classMap[label] = myClass [ both $label and $myClass are type
     * of uint64_t ]
     * */
    classifyNodes(pattern_graph, NodeSequence, classMap, parentMap, P1set,
                       S1set);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_end = MPI_Wtime();

    /********************* vf3 matching ************************/
    /*********** basic typedef **************/
    graph_type::vertex_data<uint16_t, std::allocator<uint16_t>> bfs_level_data(
        *graph);
    graph_type::vertex_data<graph_type::vertex_locator,
                            std::allocator<graph_type::vertex_locator>>
                                       bfs_parent_data(*graph);
    std::vector<std::vector<uint64_t>> resultList;
    if (mpi_rank == 0) {
      std::cout << "start matching:" << std::endl;
    }

    /*************** isomorphism matching ******************/
    MPI_Barrier(MPI_COMM_WORLD);
    double      a_start = MPI_Wtime();
    VertexClass startClass =
        classMap[pattern_graph.vertex_data[NodeSequence[0]]];  // to DO
    // debug_shell<<"[debug] start vertex:\n";
    std::vector<graph_type::vertex_locator> starting_vertices;
    int                                count = 0;
    for (auto it = graph->vertices_begin(); it != graph->vertices_end(); it++) {
      // debug_shell<<vertex_metadata[*it]<<std::endl;
      if (classMap[vertex_metadata[*it]] == startClass) {
        count++;
        // debug_shell<<graph->locator_to_label(*it)<<std::endl;
        starting_vertices.push_back(*it);
        // debug_shell<<"push succeed "<<count<<std::endl;
      }
    }
    // debug_shell<<"[debug] start matching...\n";

    // debug print
    if (mpi_rank == 0) {
      std::cout << "----------- Probabilities Map ----------------"
                << std::endl;
      for (auto it = Pf.begin(); it != Pf.end(); it++) {
        std::cout << it->first
                  << " label: " << pattern_graph.vertex_data[it->first]
                  << " value: " << it->second << std::endl;
      }
      std::cout << "-------------- pattern graph ------------------------"
                << std::endl;
      for (int v = 0; v < pattern_graph.vertex_count; v++) {
        std::cout << v << " Vertex data: " << pattern_graph.vertex_data[v]
                  << " Vertex degree: " << pattern_graph.vertex_degree[v]
                  << std::endl;
      }
      std::cout << "-------------- Node Sequence ------------------------"
                << std::endl;
      for (auto it = NodeSequence.begin(); it != NodeSequence.end(); it++) {
        std::cout << *it << " Vertex data: " << pattern_graph.vertex_data[*it]
                  << " Vertex degree: " << pattern_graph.vertex_degree[*it]
                  << std::endl;
      }
      std::cout << "-------------------- parent map ----------------------"
                << std::endl;
      for (auto it = parentMap.begin(); it != parentMap.end(); it++) {
        std::cout << it->first << " -> " << it->second << std::endl;
      }
      std::cout << "------------------ PS1 set ---------------------------"
                << std::endl;
      for (auto it = P1set.begin(); it != P1set.end(); it++) {
        std::cout << "level: " << it->first << " ";
        auto levelMap = it->second;
        for (auto it1 = levelMap.begin(); it1 != levelMap.end(); it1++) {
          std::cout << it1->first << ":{";
          auto setM = it1->second;
          for (auto it2 = setM.begin(); it2 != setM.end(); it2++) {
            std::cout << *it2 << ", ";
          }
          std::cout << "} , ";
        }
        std::cout << std::endl;
      }
    }

    vf3_pattern_matching(graph, &pattern_graph, bfs_level_data,
                              bfs_parent_data, starting_vertices, classMap,
                              NodeSequence, parentMap, vertex_metadata,
                              resultList, debug_file_name, debug_shell);
    MPI_Barrier(MPI_COMM_WORLD);
    double a_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "matching time: " << a_end - a_start << " s" << std::endl;
    }
  }
  return 0;
}
