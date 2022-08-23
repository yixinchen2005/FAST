#include <bitset>
#include <string>

#include <havoqgt/delegate_partitioned_graph.hpp>
#include <havoqgt/distributed_db.hpp>
#include <metadata/vertex_data_db.hpp>
#include <metadata/vertex_data_db_degree.hpp>
#include <prunejuice/algorithm_state_simplified.hpp>
#include <prunejuice/local_filtering.hpp>
#include <prunejuice/template_simplified.hpp>
//#include <prunejuice/local_constraint_checking.hpp>

using namespace havoqgt;

typedef delegate_partitioned_graph<distributed_db::allocator<>> graph_type;

void usage() {
  if (comm_world().rank() == 0) {
    std::cerr
        << "Usage: -i <string> -p <string> -o <string>\n"
        << " -i <string>   - input graph base filename (required)\n"
        << " -b <string>   - backup graph base filename. If set, \"input\" "
           "graph will be deleted if it exists\n"
        << " -v <string>   - vertex metadata base filename (optional, Default "
           "is degree based metadata)\n"
        << " -p <string>   - pattern base directory (required)\n"
        << " -o <string>   - output base directory (required)\n"
        << " -u <bool>     - Treat edgelist as undirected (Default is 0)\n"
        << " -h            - print help and exit\n\n";
  }
}

void parse_cmd_line(int argc, char** argv, std::string& graph_input,
                    std::string& backup_graph_input,
                    std::string& vertex_metadata_input,
                    std::string& pattern_input, bool& undirected,
                    std::string& result_output) {
  if (comm_world().rank() == 0) {
    std::cout << "CMD Line: ";
    for (int i = 0; i < argc; i++) {
      std::cout << " " << argv[i];
    }
    std::cout << std::endl;
  }

  undirected                = false;
  bool           print_help = false;
  std::bitset<2> required_input;
  required_input.reset();

  char c;
  while ((c = getopt(argc, argv, "i:b:v:p:o:u:h ")) != -1) {
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
      case 'p':
        pattern_input = optarg;
        break;
      case 'o':
        result_output = optarg;
        required_input.set(1);
        break;
      case 'u':
        undirected = atoi(optarg);
        break;
      default:
        std::cerr << "Unrecognized option: " << c << ", Ignore." << std::endl;
        print_help = true;
        break;
    }
  }

  if (print_help || !required_input.all()) {
    usage();
    exit(-1);
  }
}

int main(int argc, char** argv) {
  int mpi_rank(0), mpi_size(0);
  havoqgt::init(&argc, &argv);
  {
    CHK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    CHK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));

    if (mpi_rank == 0) {
      std::cout << "MPI Initialized with " << mpi_size << " Ranks" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string graph_input;
    std::string backup_graph_input;
    std::string vertex_metadata_input;
    std::string pattern_dir;
    bool        undirected;
    std::string result_dir;

    parse_cmd_line(argc, argv, graph_input, backup_graph_input,
                   vertex_metadata_input, pattern_dir, undirected, result_dir);
    MPI_Barrier(MPI_COMM_WORLD);

    // Loading graph
    if (mpi_rank == 0) {
      std::cout << "Loading Graph ..." << std::endl;
    }

    if (backup_graph_input.size() > 0) {
      distributed_db::transfer(backup_graph_input.c_str(), graph_input.c_str());
    }

    distributed_db ddb(db_open_read_only(), graph_input.c_str());

    auto graph = ddb.get_manager()->find<graph_type>("graph_obj").first;
    assert(graph != nullptr);

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

    typedef uint64_t Vertex;
    typedef uint64_t VertexData;
    typedef uint64_t Edge;
    // typedef uint8_t EdgeData;

    static constexpr size_t max_bit_vector_size = MAX_PATTERN_SIZE;
    // static constexpr size_t max_bit_vector_size = 32;
    typedef std::bitset<max_bit_vector_size>         BitSet;
    typedef std::unordered_map<VertexData, uint64_t> NLFMap;

    typedef graph_type::vertex_iterator vitr_type;
    typedef graph_type::vertex_locator  vloc_type;

    typedef graph_type::vertex_data<VertexData, std::allocator<VertexData>>
        VertexMetaData;

    typedef prunejuice::vertex_state<BitSet, NLFMap> VertexState;
    typedef std::unordered_map<Vertex, VertexState>  VertexStateMap;

    typedef uint8_t Boolean;
    typedef graph_type::vertex_data<Boolean, std::allocator<Boolean>>
        VertexActive;

    typedef std::unordered_map<Vertex, uint8_t> VertexUint8Map;
    typedef graph_type::vertex_data<VertexUint8Map,
                                    std::allocator<VertexUint8Map>>
        VertexUint8MapCollection;

    typedef uint16_t TemplateVertexType;
    typedef graph_type::vertex_data<TemplateVertexType,
                                    std::allocator<TemplateVertexType>>
        TemplateVertex;

    if (mpi_rank == 0) {
      std::cout << "Pattern Matching ... " << std::endl;
    }

    /* Declare vertex containers */
    VertexStateMap           vertex_state_map;
    VertexMetaData           vertex_metadata(*graph);
    VertexActive             vertex_active(*graph);
    VertexUint8MapCollection vertex_active_edges_map(*graph);
    TemplateVertex           template_candidates(*graph);

    if (mpi_rank == 0) {
      std::cout << "Pattern Matching | Allocated Vertex and Edge Containers"
                << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (vertex_metadata_input.size() > 0) {
      vertex_data_db_nostdfs<graph_type, VertexMetaData, Vertex, VertexData>(
          graph, vertex_metadata, vertex_metadata_input, 10000);
    } else {
      vertex_data_db_degree<graph_type, VertexMetaData, Vertex, VertexData>(
          graph, vertex_metadata);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    bool do_output_vertex_data = true;

    /* Output vertices of the data graph */
    if (do_output_vertex_data) {
      std::string vertex_data_filename = result_dir +
                                         "/all_ranks_vertex_data/vertex_data_" +
                                         std::to_string(mpi_rank);
      std::ofstream vertex_data_file(vertex_data_filename, std::ofstream::out);

      for (vitr_type vitr = graph->vertices_begin();
           vitr != graph->vertices_end(); ++vitr) {
        vloc_type vertex = *vitr;
        vertex_data_file << mpi_rank << ", l, "
                         << graph->locator_to_label(vertex) << ", "
                         << vertex_metadata[vertex] << ", "
                         << graph->degree(vertex) << 
                         ", " << graph->outgoing_degree(vertex) << ", "
                         << graph->incoming_degree(vertex) << "\n";
      }

      for (vitr_type vitr = graph->delegate_vertices_begin();
           vitr != graph->delegate_vertices_end(); ++vitr) {
        vloc_type   vertex = *vitr;
        std::string role =
            (vertex.is_delegate() && (graph->master(vertex) == mpi_rank))
                ? ", c, "
                : ", d, ";
        vertex_data_file << mpi_rank << role << graph->locator_to_label(vertex)
                         << ", " << vertex_metadata[vertex] << ", "
                         << graph->degree(vertex) << "\n";
      }

      vertex_data_file.close();
    }  // end of if

    if (mpi_rank == 0) {
      std::cout << "Setting up Pattern ... " << std::endl;
    }

    typedef prunejuice::pattern_graph_csr<Vertex, Edge, VertexData>
        PatternGraph;

    std::string  pattern_input_basename = pattern_dir + "/pattern";
    PatternGraph pattern_graph(pattern_input_basename + "_edge",
                               pattern_input_basename + "_vertex_data",
                               !undirected, false);

    if (mpi_rank == 0) {
      std::cout << "Pattern Matching | Searching Pattern " << std::endl;
      std::cout << "Directed: " << pattern_graph.directed << std::endl;
      for (auto v = 0; v < pattern_graph.vertex_count; v++) {
        std::cout << v << ": off-set " << pattern_graph.vertices[v]
                  << " vertex_data " << pattern_graph.vertex_data[v]
                  << " vertex_degree " << pattern_graph.vertex_degree[v]
                  << " vertex_outgoing_degree " << pattern_graph.vertex_degree[v]
                  << " vertex_incoming_degree " << pattern_graph.vertex_incoming_degree[v]
                  << std::endl;
        std::cout << " neighbors: ";
        for (auto e = pattern_graph.vertices[v];
             e < pattern_graph.vertices[v + 1]; e++) {
          auto v_nbr = pattern_graph.edges[e];
          std::cout << v_nbr << ", ";
        }
        std::cout << std::endl;
      }
    }

    /* Initialize vertex containers */
    vertex_state_map.clear();
    vertex_active.reset(true);
    vertex_active_edges_map.clear();

    /* Initialize application parameters */
    bool global_init_step    = true;
    bool global_not_finished = false;

    /* Initialized result file handlers */
    std::string   superstep_result_filename = result_dir + "/result_superstep";
    std::ofstream superstep_result_file(superstep_result_filename,
                                        std::ofstream::out);
    std::string   message_count_result_filename =
        result_dir + "/all_ranks_messages/messages_" + std::to_string(mpi_rank);
    std::ofstream message_count_result_file(message_count_result_filename,
                                            std::ofstream::out);
    std::string   active_vertices_count_result_filename =
        result_dir + "/all_ranks_active_vertices_count/active_vertices_" +
        std::to_string(mpi_rank);
    std::ofstream active_vertices_count_result_file(
        active_vertices_count_result_filename, std::ofstream::out);
    std::string active_vertices_result_filename =
        result_dir + "/all_ranks_active_vertices/active_vertices_" +
        std::to_string(mpi_rank);
    std::ofstream active_vertices_result_file(active_vertices_result_filename,
                                              std::ofstream::out);
    std::string   active_edges_count_result_filename =
        result_dir + "/all_ranks_active_edges_count/active_edges_" +
        std::to_string(mpi_rank);
    std::ofstream active_edges_count_result_file(
        active_edges_count_result_filename, std::ofstream::out);
    std::string active_edges_result_filename =
        result_dir + "/all_ranks_active_edges/active_edges_" +
        std::to_string(mpi_rank);
    std::ofstream active_edges_result_file(active_edges_result_filename,
                                           std::ofstream::out);

    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) {
      std::cout << "Running Constraint Checking ..." << std::endl;
    }

    global_not_finished = false;

    double label_propagation_time_start = MPI_Wtime();

    prunejuice::label_propagation_label_degree_filtering_bsp<
        Vertex, VertexData, VertexMetaData, BitSet, NLFMap, graph_type,
        PatternGraph, VertexStateMap, VertexActive, VertexUint8MapCollection,
        TemplateVertex>(
        graph, vertex_metadata, pattern_graph, vertex_state_map, vertex_active,
        vertex_active_edges_map, template_candidates, global_init_step,
        global_not_finished, superstep_result_file, message_count_result_file,
        active_vertices_count_result_file, active_edges_count_result_file); 

    MPI_Barrier(MPI_COMM_WORLD);
    double label_propagation_time_end = MPI_Wtime();
    if (mpi_rank == 0) {
      std::cout << "Pattern Matching Time | Local Constraint Checking : "
                << label_propagation_time_end - label_propagation_time_start
                << std::endl;
    }

    if (global_init_step) {  // Important
      global_init_step = false;
    }

    // global termination detection
    global_not_finished = havoqgt::mpi_all_reduce(
        global_not_finished, std::greater<uint8_t>(), MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0) {
      std::cout << "Pattern Matching | Global Finished Status : ";
      if (global_not_finished) {
        std::cout << "Continue" << std::endl;
      } else {
        std::cout << "Stop" << std::endl;
      }
    }

    // Important : This may slow things down -only for presenting results

    bool present_results = true;

    if (present_results) {
      for (auto& v : vertex_state_map) {
        auto v_locator = graph->label_to_locator(v.first);
        if (v_locator.is_delegate() && (graph->master(v_locator) == mpi_rank)) {
          BitSet v_template_candidates(v.second.template_candidates);
          active_vertices_result_file << mpi_rank << ", " << v.first << ", "
                                      << vertex_metadata[v_locator] << ", "
                                      << v_template_candidates << std::endl;

          // edges
          for (auto& n : vertex_active_edges_map[v_locator]) {
            active_edges_result_file << mpi_rank << ", " << v.first << ", "
                                     << n.first << std::endl;
          }

        } else if (!v_locator.is_delegate()) {
          BitSet v_template_candidates(template_candidates[v_locator]);
          active_vertices_result_file << mpi_rank << ", " << v.first << ", "
                                      << vertex_metadata[v_locator] << ", "
                                      << v_template_candidates << std::endl;

          // edges
          for (auto& n : vertex_active_edges_map[v_locator]) {
            active_edges_result_file << mpi_rank << ", " << v.first << ", "
                                     << n.first << std::endl;
          }
        }
      }
    }

    superstep_result_file.close();
    message_count_result_file.close();
    active_vertices_count_result_file.close();
    active_edges_count_result_file.close();
    active_vertices_result_file.close();
    active_edges_result_file.close();
  }  // end of havoqgt block
  return 0;
}