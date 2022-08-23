#pragma once
//#include <algorithm>
//#include <fstream>
//#include <set>
//#include <sstream>
//#include <string>
//#include <unordered_map>
//#include <vector>

//#include <havoqgt/breadth_first_search.hpp>
//#include <havoqgt/cache_utilities.hpp>
//#include <havoqgt/delegate_partitioned_graph.hpp>
//#include <havoqgt/gen_preferential_attachment_edge_list.hpp>

#include <boost/bimap.hpp>
//#include <boost/interprocess/managed_heap_memory.hpp>

#define NULL_NODE 0xffffffffffffffff

using namespace havoqgt;

namespace prunejuice
{
  typedef havoqgt::distributed_db::segment_manager_type segment_manager_t;
  typedef havoqgt::delegate_partitioned_graph<
      typename segment_manager_t::template allocator<void>::type>
      graph_type;
  typedef graph_type::vertex_locator Vertex_locator;
  typedef uint64_t VertexData;
  typedef uint64_t Vertex;
  typedef uint64_t VertexClass;
  typedef int Count; // to DO make it support the type $(unsigned long long),
                     // potential risks of overflowing
  typedef boost::bimap<Vertex, Vertex_locator> bm_type;

  template <typename PatternGraph, typename DataGraph, typename ClassMap,
            typename ParentMap, typename VertexMetadata, typename NodeSequence>
  class VF3_state
  {
  public:
    boost::bimap<Vertex, Vertex_locator> content;

    std::unordered_map<VertexClass, std::set<Vertex>> P1_ci;
    std::unordered_map<VertexClass, std::set<Vertex_locator>> P2_ci;
    int size;

    VF3_state() : size(0) {}

    VF3_state(std::vector<Vertex_locator> vL, DataGraph *dataGraph,
              PatternGraph *patternGraph, ClassMap &classMap,
              VertexMetadata &dataGraphLabel, NodeSequence &NG1)
    {
      int index = 0;
      size = 0;
      for (auto it = vL.begin(); it != vL.end(); it++)
      {
        push_pair(dataGraph, patternGraph, classMap, dataGraphLabel, NG1[index],
                  *it);
        index++;
      }
    }

    void push_pair(DataGraph *dataGraph, PatternGraph *patternGraph,
                   ClassMap &classMap, VertexMetadata &dataGraphLabel, Vertex u,
                   Vertex_locator v)
    {
      gen_P1_ci(patternGraph, classMap, u, v);
      gen_P2_ci(dataGraph, classMap, u, v, dataGraphLabel);
      // sequence cannot be changed
      content.insert(bm_type::value_type(u, v));
      size++;
    }
    bool is_feasible(DataGraph *dataGraph, PatternGraph *patternGraph,
                     ClassMap &classMap, Vertex u, Vertex_locator v,
                     bool wh_restrict)
    {
      /******** Function Fs *********/
      std::set<Vertex> P1_u;
      std::set<Vertex_locator> P2_v;
      P1(patternGraph, u, P1_u);
      P2(dataGraph, v, P2_v);
      if (!Fci(P1_u, P2_v, false))
      {
        return false;
      }
      if (wh_restrict && !Fci(P1_u, P2_v, true))
      {
        return false;
      }
      // function Fal1 and Fal2
      std::set<VertexClass> classSet;
      for (auto it = classMap.begin(); it != classMap.end(); it++)
      {
        classSet.insert(it->second);
      }
      for (auto it = classSet.begin(); it != classSet.end(); it++)
      {
        if (!Fla12i(*it, P1_u, P2_v))
        {
          return false;
        }
      }
      return true;
    }
    // just for test
    void generate_candidate(DataGraph *dataGraph, PatternGraph *patternGraph,
                            ParentMap &parentMap, Vertex u,
                            VertexMetadata &dataGraphLabel,
                            std::set<uint64_t> &candidateResult)
    {
      auto tmp = parentMap.find(u);
      if (tmp == parentMap.end())
      { // has no parent
        // std::cout<<"No parents!"<<std::endl;//to
        // DO this may cause some solutions missed
      }
      else
      {
        Vertex u_parent = tmp->second;
        Vertex_locator v_parent = (content.left.find(u_parent))->second;
        for (auto ite = dataGraph->edges_begin(v_parent);
             ite != dataGraph->edges_end(v_parent); ite++)
        {
          auto v_neighbor = ite.target();
          if (dataGraphLabel[v_neighbor] == patternGraph->vertex_data[u])
          {
            candidateResult.insert(dataGraph->locator_to_label(v_neighbor));
          }
        }
      }
    }

    void print_myself(DataGraph *dataGraph, PatternGraph *patternGraph)
    {
      std::cout << "solution: ";
      for (auto it = content.left.begin(); it != content.left.end(); it++)
      {
        std::cout << it->first << " -> "
                  << dataGraph->locator_to_label(it->second) << " , ";
      }
      std::cout << std::endl;
    }

    VF3_state &operator=(const VF3_state &other)
    {
      content = other.content;
      P1_ci = other.P1_ci;
      P2_ci = other.P2_ci;
      return *(this);
    }

  private:
    // data graph
    void P2(DataGraph *dataGraph, Vertex_locator v,
            std::set<Vertex_locator> &result)
    {
      for (auto ite = dataGraph->edges_begin(v); ite != dataGraph->edges_end(v);
           ite++)
      {
        Vertex_locator v_ = ite.target();
        result.insert(v_);
      }
    }
    // pattern graph
    void P1(PatternGraph *patternGraph, Vertex u, std::set<Vertex> &result)
    {
      for (auto e = patternGraph->vertices[u]; e < patternGraph->vertices[u + 1];
           e++)
      {
        Vertex u_ = patternGraph->edges[e];
        result.insert(u_);
      }
    }
    void gen_P1_ci(PatternGraph *patternGraph, ClassMap &classMap, Vertex u,
                   Vertex_locator v)
    {
      for (auto e = patternGraph->vertices[u]; e < patternGraph->vertices[u + 1];
           e++)
      {
        Vertex u_ = patternGraph->edges[e];
        VertexClass u_class = classMap[patternGraph->vertex_data[u_]];
        if (content.left.find(u_) == content.left.end())
        {
          if (P1_ci.find(u_class) == P1_ci.end())
          {
            std::set<Vertex> tmp;
            tmp.insert(u_);
            P1_ci[u_class] = tmp;
          }
          else
          {
            P1_ci[u_class].insert(u_);
          }
        }
      }
      // remove the vertex u if it is already in the $P1_ci
      VertexClass p_class = classMap[patternGraph->vertex_data[u]];
      auto it = P1_ci.find(p_class);
      if (it != P1_ci.end())
      {
        auto it_ = P1_ci[p_class].find(u);
        if (it_ != P1_ci[p_class].end())
        {
          P1_ci[p_class].erase(it_);
        }
      }
    }
    void gen_P2_ci(DataGraph *dataGraph, ClassMap &classMap, Vertex u,
                   Vertex_locator v, VertexMetadata &dataGraphLabel)
    {
      for (auto ite = dataGraph->edges_begin(v); ite != dataGraph->edges_end(v);
           ite++)
      {
        Vertex_locator v_ = ite.target();
        VertexClass v_class = classMap[dataGraphLabel[v_]];
        if (content.right.find(v_) == content.right.end())
        {
          if (P2_ci.find(v_class) == P2_ci.end())
          {
            std::set<Vertex_locator> tmp;
            tmp.insert(v_);
            P2_ci[v_class] = tmp;
          }
          else
          {
            P2_ci[v_class].insert(v_);
          }
        }
      }
      // remove the vertex v if it is already in the $P2_ci
      VertexClass d_class = dataGraphLabel[v];
      auto it = P2_ci.find(d_class);
      if (it != P2_ci.end())
      {
        auto it_ = P2_ci[d_class].find(v);
        if (it_ != P2_ci[d_class].end())
        {
          P2_ci[d_class].erase(it_);
        }
      }
    }

    bool Fci(std::set<Vertex> &P1_u, std::set<Vertex_locator> &P2_v,
             bool reverse)
    {
      for (auto it = content.begin(); it != content.end(); it++)
      {
        if (!reverse)
        {
          Vertex u = it->left;
          if (P1_u.find(u) != P1_u.end())
          {
            Vertex_locator v = it->right;
            if (P2_v.find(v) == P2_v.end())
            {
              std::cout << " ooo1ooo ";
              return false;
            }
          }
        }
        else
        {
          Vertex_locator v = it->right;
          if (P2_v.find(v) != P2_v.end())
          {
            Vertex u = it->left;
            if (P1_u.find(u) == P1_u.end())
            {
              std::cout << " ooo2ooo ";
              return false;
            }
          }
        }
      }
      return true;
    }
    bool Fla12i(VertexClass vertexClass, std::set<Vertex> &P1_u,
                std::set<Vertex_locator> &P2_v)
    {
      int count_left = 0, count_right = 0;
      int count_left_v = 0, count_right_v = 0;
      for (auto it = P1_u.begin(); it != P1_u.end(); it++)
      {
        auto P1set_ci = P1_ci[vertexClass];
        if (P1set_ci.find(*it) != P1set_ci.end())
        {
          count_left++;
        }
        else
        {
          count_left_v++;
        }
      }
      for (auto it = P2_v.begin(); it != P2_v.end(); it++)
      {
        auto P2set_ci = P2_ci[vertexClass];
        if (P2set_ci.find(*it) != P2set_ci.end())
        {
          count_right++;
        }
        else
        {
          count_right_v++;
        }
      }
      if (count_left > count_right)
      {
        std::cout << " ooo3ooo ";
      }
      if (count_left_v > count_left_v)
      {
        std::cout << " ooo3ooo ";
      }
      return (count_left <= count_right) && (count_left_v <= count_right_v);
    }
  };

  template <typename Graph, typename State>
  class vf3_visitor
  {
  public:
    typedef typename Graph::vertex_locator vertex_locator;
// contruction function
#pragma GCC diagnostic ignored \
    "-Woverflow" // NOTE: is there a better way to clean these overflows?
    vf3_visitor() : m_level(std::numeric_limits<uint64_t>::max())
    {
    }
#pragma GCC diagnostic pop
    vf3_visitor(Vertex_locator _vertex, uint64_t _level, Vertex_locator _parent,
                std::vector<vertex_locator> partial_match,
                std::vector<uint64_t> partial_label, int msg)
        : vertex(_vertex), m_parent(_parent), m_level(_level), message(msg)
    {
      int i = 0;
      for (auto it = partial_match.begin(); it != partial_match.end(); it++)
      {
        set_match(i, *it);
        set_label(i, partial_label[i]);
        i++;
      }
    }
    vf3_visitor(Vertex_locator _vertex)
        : vertex(_vertex), m_parent(_vertex), m_level(1), message(0) {}

    template <typename AlgData>
    bool pre_visit(AlgData &alg_data) const
    {
      return true;
    }

    template <typename VisitorQueueHandle, typename AlgData>
    bool init_visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
    {
      return visit(g, vis_queue, alg_data);
    }

    vertex_locator get_match(int i)
    {
      switch (i)
      {
      case 0:
        return v0;
      case 1:
        return v1;
      case 2:
        return v2;
      case 3:
        return v3;
      case 4:
        return v4;
      case 5:
        return v5;
      case 6:
        return v6;
      case 7:
        return v7;
      }
    }

    void set_match(int i, vertex_locator v)
    {
      switch (i)
      {
      case 0:
        v0 = v;
        return;
      case 1:
        v1 = v;
        return;
      case 2:
        v2 = v;
        return;
      case 3:
        v3 = v;
        return;
      case 4:
        v4 = v;
        return;
      case 5:
        v5 = v;
        return;
      case 6:
        v6 = v;
        return;
      case 7:
        v7 = v;
        return;
      }
    }

    uint64_t get_label(int i)
    {
      switch (i)
      {
      case 0:
        return l0;
      case 1:
        return l1;
      case 2:
        return l2;
      case 3:
        return l3;
      case 4:
        return l4;
      case 5:
        return l5;
      case 6:
        return l6;
      case 7:
        return l7;
      }
    }

    void set_label(int i, uint64_t l)
    {
      switch (i)
      {
      case 0:
        l0 = l;
        return;
      case 1:
        l1 = l;
        return;
      case 2:
        l2 = l;
        return;
      case 3:
        l3 = l;
        return;
      case 4:
        l4 = l;
        return;
      case 5:
        l5 = l;
        return;
      case 6:
        l6 = l;
        return;
      case 7:
        l7 = l;
        return;
      }
    }

    void update()
    {
      partial_match.clear();
      partial_label.clear();
      for (int i = 0; i < level(); i++)
      {
        partial_match.push_back(get_match(i));
        partial_label.push_back(get_label(i));
      }
    }

    template <typename AlgData>
    bool isFeasible(Graph &g, Vertex current_un, vertex_locator candidate_vn,
                    AlgData &alg_data)
    {
      if (partial_match.size() == 1)
      {
        return true;
      }
      auto patternGraph = std::get<2>(alg_data);
      auto NG1 = std::get<4>(alg_data);
      auto vertexMetadata = std::get<6>(alg_data);
      if (patternGraph->vertex_data[current_un] != vertexMetadata[candidate_vn])
      {
        return false;
      }

      for (int i = 0; i < partial_match.size() - 1; i++)
      {
        if (g.locator_to_label(candidate_vn) == get_label(i))
        {
          return false;
        }
      }

      for (auto e = patternGraph->vertices[current_un];
           e < patternGraph->vertices[current_un + 1]; e++)
      {
        auto u_neighbor = patternGraph->edges[e];

        int index = 0;
        for (; index < partial_match.size(); index++)
        {
          if (NG1[index] == u_neighbor)
          {
            if (index < partial_match.size())
            {
              auto vn_label = get_label(index); // partial_match[index];
              bool found = false;
              for (auto ite = g.edges_begin(vertex); ite != g.edges_end(vertex);
                   ite++)
              {
                auto v_neighbor = g.locator_to_label(ite.target());
                if (vn_label == v_neighbor)
                {
                  found = true;
                  break;
                }
              }
              if (found == false)
              {
                return false;
              }
            }
            break;
          }
        }
      }
      return true;
    }

    template <typename AlgData>
    vertex_locator getVn(Vertex un,
                         AlgData &alg_data /*, bool debug, Graph& g*/)
    {
      auto NG1 = std::get<4>(alg_data);
      int index = 0;
      for (int index = 0; index < NG1.size(); index++)
      {
        if (NG1[index] == un)
        {
          if (index >= partial_match.size())
          {
            std::cout << "Error: index overflows" << std::endl;
          }
          return get_match(index);
        }
      }
      std::cout << "Error: not corrent" << std::endl;
      return get_match(index); // partial_match[index];
    }

    template <typename VisitorQueueHandle, typename AlgData>
    bool visit(Graph &g, VisitorQueueHandle vis_queue, AlgData &alg_data)
    {
      auto NG1 = std::get<4>(alg_data);
      auto vertexMetadata = std::get<6>(alg_data);
      auto parentMap = std::get<5>(alg_data);
      auto patternGraph = std::get<2>(alg_data);

      // debug
      int mpi_rank = havoqgt::comm_world().rank();

      if (message == 0)
      {
        // partial_match.push_back(vertex);
        set_match(level() - 1, vertex);
        set_label(level() - 1, g.locator_to_label(vertex));
        update();
        Vertex current_un = NG1[partial_match.size() - 1];
        if (isFeasible(g, current_un, vertex, alg_data))
        {
          if (partial_match.size() == NG1.size())
          {
            std::string filePath = std::get<8>(alg_data);
            std::ofstream outP(filePath, std::ofstream::app);
            for (int i = 0; i < NG1.size(); i++)
            {
              outP << get_label(i) << " ";
            }
            outP << std::endl;
            outP.close();
            return true;
          }
          Vertex next_un = NG1[partial_match.size()];
          if (parentMap.find(next_un) == parentMap.end())
          {
            for (auto it = g.vertices_begin(); it != g.vertices_end(); it++)
            {
              vf3_visitor new_visitor(*it, level() + 1, vertex, partial_match,
                                      partial_label, 0);
              vis_queue->queue_visitor(new_visitor);
            }
          }
          else
          {
            Vertex parent_un = parentMap[next_un];
            vertex_locator parent_vn = getVn(parent_un, alg_data);
            vf3_visitor new_visitor(parent_vn, level(), vertex, partial_match,
                                    partial_label, 1);
            vis_queue->queue_visitor(new_visitor);
          }
        }
      }
      else
      {
        update();
        for (auto ite = g.edges_begin(vertex); ite != g.edges_end(vertex);
             ite++)
        {
          auto v_neighbor = ite.target();
          vf3_visitor new_visitor(v_neighbor, level() + 1, vertex, partial_match,
                                  partial_label, 0);
          vis_queue->queue_visitor(new_visitor);
        }
      }
      return false;
    }

    uint64_t level() const { return m_level; }
    Vertex_locator parent() const { return m_parent; }

    friend inline bool operator>(const vf3_visitor &v1, const vf3_visitor &v2)
    {
      if (v1.level() > v2.level())
      {
        return false;
      }
      return true;
    }

    vertex_locator vertex;
    vertex_locator m_parent;
    vertex_locator v0, v1, v2, v3, v4, v5, v6, v7;
    uint64_t l0, l1, l2, l3, l4, l5, l6, l7;
    uint64_t m_level : 16;
    std::vector<vertex_locator> partial_match;
    std::vector<uint64_t> partial_label;
    int message; // 0 generate_candidate | 1 find neighbors
  };

  // starting the asynchronous algorithm
  template <typename TGraph, typename PatternGraph, typename LevelData,
            typename ParentData, typename StartVerticesList, typename ClassMap,
            typename NodeSequence, typename ParentMap, typename VertexMetadata,
            typename ResultList>
  void vf3_pattern_matching(TGraph *g, PatternGraph *patternGraph,
                            LevelData &level_data, ParentData &parent_data,
                            StartVerticesList &startVertices, ClassMap &classMap,
                            NodeSequence &nodeSequence, ParentMap &parentMap,
                            VertexMetadata &vertexMetadata, ResultList &result,
                            std::string &outputFilePath,
                            std::ofstream &debug_shell)
  {
    typedef VF3_state<PatternGraph, TGraph, ClassMap, ParentMap, VertexMetadata,
                      NodeSequence>
        State;
    typedef vf3_visitor<TGraph, State> visitor_type;
    std::vector<std::string> resultLL;
    int count = 0;
    auto alg_data = std::forward_as_tuple(level_data, parent_data, patternGraph,
                                          classMap, nodeSequence, parentMap,
                                          vertexMetadata, result, outputFilePath);
    auto vq = create_visitor_queue<visitor_type,
                                   havoqgt::detail::visitor_priority_queue>(
        g, alg_data);
    vq.init_visitor_traversal(startVertices);
    // debug_shell<<count<<std::endl;
    for (auto it = resultLL.begin(); it != resultLL.end(); it++)
    {
      debug_shell << *it << std::endl;
    }
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
    for (auto vitr = dataGraph->vertices_begin();
         vitr != dataGraph->vertices_end(); vitr++)
    {
      // process the degree
      if ((uint32_t)mpi_rank ==
          (*vitr).owner())
      { // if the vertex is in the process
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

    // process the unordered_map into a struct list, in order to use the MPI send
    // outgoing degree, num
    std::vector<std::pair<Count, Count>> sendBufOutDegree;
    for (auto it = outDegreeMap.begin(); it != outDegreeMap.end(); it++)
    {
      sendBufOutDegree.push_back(std::pair<Count, Count>(it->first, it->second));
    }
    // start to send the partial map
    std::vector<std::pair<Count, Count>> recvBufOutDegree;
    mpi_all_gather(
        sendBufOutDegree, recvBufOutDegree,
        MPI_COMM_WORLD); // every process will have the reduced results
    // agregate the received results
    outDegreeMap.clear();
    for (auto it = recvBufOutDegree.begin(); it != recvBufOutDegree.end(); it++)
    {
      if (outDegreeMap.find(it->first) == outDegreeMap.end())
      {
        outDegreeMap[it->first] = it->second;
      }
      else
      {
        outDegreeMap[it->first] += it->second;
      }
      total_vertex_count += it->second;
    }

    // process the unordered_map into a struct list, in order to use the MPI send
    // incoming degree, num
    std::vector<std::pair<Count, Count>> sendBufInDegree;
    for (auto it = inDegreeMap.begin(); it != inDegreeMap.end(); it++)
    {
      sendBufInDegree.push_back(std::pair<Count, Count>(it->first, it->second));
    }
    // start to send the partial map
    std::vector<std::pair<Count, Count>> recvBufInDegree;
    mpi_all_gather(
        sendBufInDegree, recvBufInDegree,
        MPI_COMM_WORLD); // every process will have the reduced results
    // agregate the received results
    inDegreeMap.clear();
    for (auto it = recvBufInDegree.begin(); it != recvBufInDegree.end(); it++)
    {
      if (inDegreeMap.find(it->first) == inDegreeMap.end())
      {
        inDegreeMap[it->first] = it->second;
      }
      else
      {
        inDegreeMap[it->first] += it->second;
      }
    }

    // same procedure to process the label data
    std::vector<std::pair<Count, Count>>
        sendBufLabel; // label(this type is supposed to $VertexDdata), num
    for (auto it = labelMap.begin(); it != labelMap.end(); it++)
    {
      sendBufLabel.push_back(std::pair<VertexData, Count>(it->first, it->second));
    }
    std::vector<std::pair<Count, Count>> recvBufLabel;
    mpi_all_gather(sendBufLabel, recvBufLabel, MPI_COMM_WORLD);
    labelMap.clear();
    for (auto it = recvBufLabel.begin(); it != recvBufLabel.end(); it++)
    {
      if (labelMap.find(it->first) == labelMap.end())
      {
        labelMap[it->first] = it->second;
      }
      else
      {
        labelMap[it->first] += it->second;
      }
    }

    // debug to print the result map
    if (mpi_rank == 0)
    {
      std::cout << "--------- total vertex count " << total_vertex_count
                << " ---------------" << std::endl;
      for (auto it = labelMap.begin(); it != labelMap.end(); it++)
      {
        std::cout << it->first << " -> " << it->second << std::endl;
      }
      std::cout << "-------------- delimiter ------------" << std::endl;
    }
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