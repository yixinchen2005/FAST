#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "omp.h"

#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "../constraint/template_constraint.hpp"
#include "../constraint/template_graph.hpp"
#include "../constraint/combination.hpp"
#include "../constraint/graph_algorithm.hpp"
#include "../constraint/edge_dfs.hpp"
#include "../constraint/file_utilities.hpp"
#include "../constraint/utilities.hpp"

#include "config.hpp"

// #define ENABLE_CONSTRAINT
#define ENABLE_ENUMERATION

template <typename T> 
bool has_common_element(T& a, T& b) {
  bool found_common_element = false;
  for (auto& i : b) {
    auto find_key = a.find(i.first);
    if (find_key != a.end()) {
      found_common_element = true;
      break; 
    } 
  }  
  return found_common_element;
}

////////////////////////////////////////////////////////////////////////////////
 
template <typename Vertex, typename Edge, typename EdgeListMap, 
  typename EdgeSet>
bool is_subset_edgeset(EdgeListMap& a, EdgeSet& b) {
  bool is_subset = true;
  for(auto it = b.begin(); it != b.end(); ++it) {
    auto find_it = a.find(it->first);  
    if (find_it == a.end()) {
      is_subset = false;
      break;   
    }
  } // for
  return is_subset;  
}  

////////////////////////////////////////////////////////////////////////////////

constexpr size_t VERTEX_COUNT = 64;
typedef std::bitset<VERTEX_COUNT> VertexBitSet;

template <typename Vertex, typename VertexBitSet>
void update_vertex_bitset(Vertex s, Vertex t, VertexBitSet& vertex_bitset) {
  vertex_bitset.set(static_cast<size_t>(s));
  vertex_bitset.set(static_cast<size_t>(t));       
}

template <typename Vertex, typename Edge>
Edge edge_hash_vertex_bitset(Vertex s, Vertex t) {
  VertexBitSet vertex_bitset; 
  vertex_bitset.set(static_cast<size_t>(s));
  vertex_bitset.set(static_cast<size_t>(t));
  return static_cast<Edge>(vertex_bitset.to_ullong());   
}

////////////////////////////////////////////////////////////////////////////////

/**
 * All template constraints
 */
template <typename Vertex, typename Edge, typename EdgeSet, 
  typename TemplateConstraint, typename GraphNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesAllPaths, typename EdgeListTuple,
  typename EdgeListTupleVector, typename EdgeListTupleVectors,
  typename TemplateConstraintVector>
void populate_template_constraints(
  GraphNonLocalPropertiesUnique& graph_cycles_unique, 
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths,
  GraphNonLocalPropertiesUnique& graph_tds_cycles_unique,
  EdgeListTupleVectors& graph_tds_cycles_all_paths,
  EdgeListTuple& graph_enumeration_unique,   
  TemplateConstraintVector& template_constraints) {
  
  for (auto it = graph_cycles_unique.begin(); 
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first, 
      (graph_cycles_all_paths.find(it->first))->second, // EdgeSetVector
      template_constraints.size(), TemplateConstraint::CYCLE));    
  } // for
  
  for (auto it = graph_tds_cycles_unique.begin();
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first,
      graph_tds_cycles_all_paths[it->first], // EdgeListTupleVector 
      template_constraints.size(), TemplateConstraint::TDS));
  } // for 

  EdgeListTupleVector graph_enumeration_path(0); // TODO: improve 
    graph_enumeration_path.push_back(graph_enumeration_unique);

  template_constraints.push_back(TemplateConstraint(0, 
    graph_enumeration_path, // EdgeListTupleVector //  
    template_constraints.size(), TemplateConstraint::ENUMERATION));  
} 

template <typename Vertex, typename Edge, typename EdgeSet, 
  typename TemplateConstraint, typename GraphNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesAllPaths, typename EdgeListTupleVectors,
  typename TemplateConstraintVector>
void populate_template_constraints(
  GraphNonLocalPropertiesUnique& graph_cycles_unique, 
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths,
  GraphNonLocalPropertiesUnique& graph_tds_cycles_unique,
  EdgeListTupleVectors& graph_tds_cycles_all_paths,   
  TemplateConstraintVector& template_constraints) {
  
  for (auto it = graph_cycles_unique.begin(); 
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first, 
      (graph_cycles_all_paths.find(it->first))->second, // EdgeSetVector
      template_constraints.size(), TemplateConstraint::CYCLE));    
  } // for
  
  for (auto it = graph_tds_cycles_unique.begin();
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first,
      graph_tds_cycles_all_paths[it->first], // EdgeListTupleVector 
      template_constraints.size(), TemplateConstraint::TDS));
  } // for  
} 

template <typename Vertex, typename Edge, typename EdgeSet,
  typename TemplateConstraint, typename GraphNonLocalPropertiesUnique,
  typename GraphNonLocalPropertiesAllPaths,
  typename TemplateConstraintVector>
void populate_cycle_constraints(
  GraphNonLocalPropertiesUnique& graph_cycles_unique,
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths,
  TemplateConstraintVector& template_constraints) {

  for (auto it = graph_cycles_unique.begin();
    it != graph_cycles_unique.end(); ++it) {
    template_constraints.push_back(TemplateConstraint(it->second, it->first,
      (graph_cycles_all_paths.find(it->first))->second, // EdgeSetVector
      template_constraints.size(), TemplateConstraint::CYCLE));
  } // for
}

////////////////////////////////////////////////////////////////////////////////

template <typename Vertex, typename VertexSet, typename EdgeSet>
VertexSet edgeset_to_vertexset(EdgeSet& edgeset) {
  VertexSet vertexset(0);
  for (auto e : edgeset) {
    auto find_s = vertexset.find(std::get<0>(e.second)); 
    if (find_s == vertexset.end()) {
      vertexset.insert(std::get<0>(e.second));
    }
    
    auto find_t = vertexset.find(std::get<1>(e.second));  
    if (find_t == vertexset.end()) {
      vertexset.insert(std::get<1>(e.second));
    }
  } 
  return vertexset; 
}

/**
 * Generate TDS cycle constraints
 */
template <typename Vertex, typename Edge, typename EdgeSet,
  typename GraphNonLocalPropertiesUnique, typename  EdgeListTuple, 
  typename EdgeListTupleVector, typename EdgeListTupleVectors>
void generate_tds_cycle_constraints
  (GraphNonLocalPropertiesUnique& graph_tds_cycles_unique, 
  EdgeListTupleVectors& graph_tds_cycles_all_paths, bool do_print = false) {
 
  typedef std::unordered_set<Vertex> VertexSet;

  GraphNonLocalPropertiesUnique new_graph_tds_cycles_unique(0);

  size_t edgeset_hash = 0;

  for (auto it = graph_tds_cycles_unique.begin(); 
    it != graph_tds_cycles_unique.end(); ++it) {
    // it->second is of type EdgeSet
    VertexSet vertexset = edgeset_to_vertexset<Vertex, VertexSet>(it->second);   

    EdgeListTupleVector edgeset_all_paths(0); 

    for (auto v : vertexset) {
      edgeset_all_paths.push_back(prunejuice::pattern::graphalgorithm::edge_dfs
        <Vertex, Edge, EdgeListTuple, EdgeSet>(it->second, v));

      // Test 
      assert(edgeset_all_paths.size() > 0);   
    } // for 
    
    graph_tds_cycles_all_paths.push_back(edgeset_all_paths);

    // update graph_tds_cycles_unique to have new edgeset_hash
    // graph_tds_cycles_all_paths[0] contains the paths for edgeset_hash == 0 
    new_graph_tds_cycles_unique.insert({edgeset_hash, it->second});
    edgeset_hash++;
  } // for

  graph_tds_cycles_unique = new_graph_tds_cycles_unique; // Important:

} 

/**
 * Generate TDS cycles from all the unique cycles
 */ 
template <typename TDSEdgeSet, typename EdgeSetIterator>
void generate_tds_cycles_recursive(TDSEdgeSet& tds_edge_set,
  EdgeSetIterator it, EdgeSetIterator it_end) {
  if (it == it_end) {
    return;
  } else {
    // add edge
    auto find_edge_hash = tds_edge_set.find(it->first);
    if (find_edge_hash == tds_edge_set.end()) {
      tds_edge_set.insert({it->first, it->second});  
    }
    generate_tds_cycles_recursive(tds_edge_set, ++it, it_end);
  }  
}

template <typename Veretx, typename Edge, 
  typename GraphNonLocalPropertiesUnique>
void generate_tds_cycles(GraphNonLocalPropertiesUnique graph_cycles_unique, // copy
  GraphNonLocalPropertiesUnique& tds_cycles) {
  
  for (auto it = graph_cycles_unique.begin(); it != graph_cycles_unique.end();) {
   
    // add to tds_cycles  
    auto find_it = tds_cycles.find(it->first); 
    if (find_it == tds_cycles.end()) {
      auto insert_status = tds_cycles.insert({it->first, it->second});

      if(!insert_status.second) {
        std::cerr << "Error: failed to add an element to the map. " << std::endl;
        return;
      } 

      auto& tds_cycles_it = insert_status.first; // iterator 
      it = graph_cycles_unique.erase(it); // C++11 // remove from graph_cycles_unique      
      bool is_tds_cycle = false;

      // find intersection of cycle constraints
      for (auto it_nest = it; it_nest != graph_cycles_unique.end();) {
 
        // if intersection is found, it_nest is merged with tds_cycles_it 
        // and removed from graph_cycles_unique  
        bool found_intersection = has_common_element
          (insert_status.first->second, it_nest->second);  

        if (found_intersection) {
          // merege cycles
          generate_tds_cycles_recursive(insert_status.first->second, 
            it_nest->second.begin(), it_nest->second.end());

          if (it == it_nest) { // Important: must update it
            it_nest = graph_cycles_unique.erase(it_nest); // C++11 // remove it_nest from graph_cycles_unique
            it = it_nest;
          } else {
            it_nest = graph_cycles_unique.erase(it_nest); // C++11 // remove it_nest from graph_cycles_unique
          }

          is_tds_cycle = true;

        } else { // no intersection found
          it_nest++;
        } 
      } // for  
      // nothing was merged with this cycle, so not a TDS cycle; 
      // therefore, remove from the list of TDS cycles 
      if (!is_tds_cycle) {
        tds_cycles.erase(tds_cycles_it);          
      }          
    } else {
      std::cerr << "Error: unexpected item in the map. " << std::endl;  
    }

  } // for 
}

//------------------------------------------------------------------------------

/**
 * Not used  
 */
template <typename TDSEdgeSet, typename EdgeSetIterator>
void generate_tds_cycles_recursive_2(TDSEdgeSet& tds_edge_set, 
  EdgeSetIterator it, EdgeSetIterator it_end) {
  if (it == it_end) {
    return;
  } else {
    // add edge
    auto find_edge_hash = tds_edge_set.find(it->first);
    if (find_edge_hash == tds_edge_set.end()) { 
      tds_edge_set.insert({it->first, std::forward_as_tuple
        (std::get<0>(it->second),std::get<1>(it->second))});
    } 
    generate_tds_cycles_recursive_2(tds_edge_set, ++it, it_end);
  }
}

template <typename Veretx, typename Edge, typename GraphNonLocalPropertiesUnique, typename TDSEdgeSet>
void generate_tds_cycles(GraphNonLocalPropertiesUnique& graph_cycles_unique, TDSEdgeSet& tds_edge_set) {
  for (auto i :  graph_cycles_unique) {
    //auto it = i.begin();
    generate_tds_cycles_recursive_2(tds_edge_set, i.second.begin(), i.second.end());
  }
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Find all the unique cycles in the graph
 */
template <typename Vertex, typename Edge, 
  typename VertexNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesUnique>
void find_unique_cycles(
  VertexNonLocalPropertiesUnique& vertex_cycles_unique, 
  GraphNonLocalPropertiesUnique& graph_cycles_unique) {

  // Note: preserving edge order in a path is not required here   

  for (size_t i = 0; i < vertex_cycles_unique.size(); i++) { // vector
    for (auto& j : vertex_cycles_unique[i]) { // j is a map       
      auto find_path = graph_cycles_unique.find(j.first); 
      if (find_path == graph_cycles_unique.end()) {
        graph_cycles_unique.insert({j.first, j.second});
      } 
    } // for  
  } // for 
}

template <typename Vertex, typename Edge,  
  typename VertexNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesUnique, 
  typename GraphNonLocalPropertiesAllPaths>
void find_unique_cycles(
  VertexNonLocalPropertiesUnique& vertex_cycles_unique, 
  GraphNonLocalPropertiesUnique& graph_cycles_unique,
  GraphNonLocalPropertiesAllPaths& graph_cycles_all_paths) {

  // Note: preserving edge order in a path is not required here   

  for (size_t i = 0; i < vertex_cycles_unique.size(); i++) { // vector
    for (auto& j : vertex_cycles_unique[i]) { // j is a map
      // j.first // path hash // Note: path hash is unique for a cycle
      // j.second // path map         
      auto find_path = graph_cycles_unique.find(j.first); 
      if (find_path == graph_cycles_unique.end()) {
        graph_cycles_unique.insert({j.first, j.second});
      }

      auto find_paths = graph_cycles_all_paths.find(j.first);
      graph_cycles_all_paths[j.first].push_back(j.second);        
    } // for  
  } // for 
}

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique,
  typename OrderedPath, typename OrderedPathEdges, typename UnorderedPathEdges>
void find_cycles_recursive(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges,
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique,  
  Vertex source_vertex, Vertex v, OrderedPath walk_history, 
  OrderedPathEdges walk_history_edges, 
  UnorderedPathEdges walk_history_edges_unordered, size_t r) {

    // v is the current vertex    
    // r must be greater than 2 for it to be considered a cycle    
    // do not forward a reference to walk_history, make a new copy

    for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
      Vertex v_nbr = edges[e];
      if (vertex_degree[v_nbr] < 2) {
        //continue;
      }
      
      if (v_nbr == source_vertex && r >= 3) { 
        // tuple - vector, size_t 
        OrderedPathEdges new_walk_history_edges(walk_history_edges);
        std::get<0>(new_walk_history_edges).push_back(std::forward_as_tuple(v, v_nbr));
        //std::get<1>(new_walk_history_edges) = 
        //  edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);
        update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges));  
        // tuple - map, size_t  
        UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);

        Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

        auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
        if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
          std::get<0>(new_walk_history_edges_unordered).
            insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered));
        }

        // add path to vertex_cycles
        // vertex_cycles_unique - vector based walk history edges     
 
        auto find_path = vertex_cycles[source_vertex].find(std::get<1>(new_walk_history_edges).to_ullong());  
        if (find_path ==  vertex_cycles[source_vertex].end()) {    
          vertex_cycles[source_vertex].
            insert({std::get<1>(new_walk_history_edges).to_ullong(), std::get<0>(new_walk_history_edges)});
        }
        // vertex_cycles_unique - map based walk history edges    
        //walk_history_edges_unordered

        auto find_path_3 = vertex_cycles_unique[source_vertex].find(std::get<1>(new_walk_history_edges_unordered).to_ullong());
        if (find_path_3 ==  vertex_cycles_unique[source_vertex].end()) {  
          vertex_cycles_unique[source_vertex].
            insert({std::get<1>(new_walk_history_edges_unordered).to_ullong(), std::get<0>(new_walk_history_edges_unordered)});    
        }

        //return; 
      } else if (v_nbr == source_vertex && r < 3) {
        // path length is invalid, ignore
        //return; 
        //continue;        
      } else if (v_nbr != source_vertex && r > 0) { 
        // forward  
        auto find_v_nbr = walk_history.find(v_nbr); 
        if (find_v_nbr == walk_history.end()) {
          //walk_history.insert({v_nbr, r});           

          // vertex -     
          OrderedPath new_walk_history(walk_history);              
          new_walk_history.insert({v_nbr, r + 1}); 
 
          //--------------------------------------------------------------------
          // tuple - vector, size_t
          OrderedPathEdges new_walk_history_edges(walk_history_edges);
          std::get<0>(new_walk_history_edges).push_back(std::forward_as_tuple(v, v_nbr));
          update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges)); 

          //-------------------------------------------------------------------- 

          // tuple - map, size_t 
          UnorderedPathEdges new_walk_history_edges_unordered(walk_history_edges_unordered);
             
          Edge edge_hash = edge_hash_vertex_bitset<Vertex, Edge>(v, v_nbr);

          auto find_edge_hash = std::get<0>(new_walk_history_edges_unordered).find(edge_hash);
          if (find_edge_hash == std::get<0>(new_walk_history_edges_unordered).end()) {
            std::get<0>(new_walk_history_edges_unordered).
              insert({edge_hash, std::forward_as_tuple(v, v_nbr, r)});//r + 1)});  
            update_vertex_bitset(v, v_nbr, std::get<1>(new_walk_history_edges_unordered)); 
          }      

          //--------------------------------------------------------------------  
          // forward
          find_cycles_recursive<Vertex, Edge, VertexList, EdgeList,
            VertexNonLocalProperties, VertexNonLocalPropertiesUnique,
            OrderedPath, OrderedPathEdges, UnorderedPathEdges>
            (vertices, vertex_count, vertex_degree, edges, vertex_cycles,
            vertex_cycles_unique,  
            source_vertex, v_nbr, new_walk_history, new_walk_history_edges, 
            new_walk_history_edges_unordered,
            r + 1);

        } else { 
          //return; // repeated vertex, ignore   
          //continue;
        }
      } // if
 
    } // for   
} 

template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename VertexNonLocalProperties, 
  typename VertexNonLocalPropertiesUnique>
void find_cycles_parallel(VertexList& vertices, Vertex vertex_count, 
  VertexList& vertex_degree, EdgeList& edges, 
  VertexNonLocalProperties& vertex_cycles, 
  VertexNonLocalPropertiesUnique& vertex_cycles_unique) {

  { 
  #pragma omp parallel for
  for (Vertex v = 0; v < vertex_count; v++) {
    //Vertex v = 0; // Test   

    if (vertex_degree[v] < 2) {
      //continue;
    }
    size_t r = 1; // walk step, initialized to 1

    // map based walk history of vertices with step order
    typedef std::unordered_map<Vertex, size_t> OrderedPath;
    OrderedPath walk_history(0); 
  
    auto find_v = walk_history.find(v); 
    if (find_v == walk_history.end()) {
      walk_history.insert({v, r});  
    } else {
      std::cerr << "Error: unexpected item in the map. " << std::endl;
    }

    // vector based walk history of edges
    //typedef std::vector<std::tuple<Vertex, Vertex>> OrderedPathEdges; 
    typedef std::tuple<std::vector<std::tuple<Vertex, Vertex>>, VertexBitSet> OrderedPathEdges;
    //OrderedPathEdges walk_history_edges(0);     
 
    //walk_history_edges.push_back(std::forward_as_tuple(v, v));
    
    std::vector<std::tuple<Vertex, Vertex>> wh_vector(0);
    //wh_vector.push_back(std::forward_as_tuple(v, v)); // not required 
    VertexBitSet vertex_bitset;
    //update_vertex_bitset(v, v, vertex_bitset); // not required
    OrderedPathEdges walk_history_edges(wh_vector, vertex_bitset);

    // map based walk history of edges with step order
    typedef std::tuple< std::unordered_map< Edge, std::tuple<Vertex, Vertex, size_t> >, VertexBitSet> 
      UnorderedPathEdges;    
    std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>> wh_unordered_map(0);
    UnorderedPathEdges walk_history_edges_unordered(wh_unordered_map, vertex_bitset); 

    find_cycles_recursive<Vertex, Edge, VertexList, EdgeList, 
      VertexNonLocalProperties, VertexNonLocalPropertiesUnique, 
      OrderedPath, OrderedPathEdges, UnorderedPathEdges>
      (vertices, vertex_count, vertex_degree, edges, vertex_cycles, 
      vertex_cycles_unique, 
      v, v, walk_history, walk_history_edges, walk_history_edges_unordered, r); 
     
  } // for
  } 
}

////////////////////////////////////////////////////////////////////////////////

/**
 * Template prototype generation - up to k-edit distance
 */

template <typename Edge, typename TemplateGraph, typename O, typename S> 
void generate_edge_combinations(TemplateGraph& input_template, 
  size_t k_edit_distance, O& edge_combinations, S& edge_set_combinations) {
  //typedef typename TemplateGraph::Edge_t Edge_t; 

  // setup input
  typedef std::vector<Edge> NCollection;
  NCollection n_items(0); 

  for (auto& e : input_template.edgelist_unique) {
    auto edge_uint = e.first;
    n_items.push_back(static_cast<Edge>(edge_uint));
  }
  std::stable_sort(n_items.begin(), n_items.end(),
    [](const Edge& a, const Edge& b) -> bool {
         return a < b;
       });

  //size_t n_collection = 7; // Test
  size_t n_collection = n_items.size();
  size_t k_selection = k_edit_distance;

  assert(n_collection >= k_selection);

  //size_t nCk = static_cast<size_t>(boost::math::binomial_coefficient<double>
  //  (static_cast<size_t>(n_collection.size()), k_selection));
  size_t nCk = static_cast<size_t>(boost::math::binomial_coefficient<double>
    (n_collection, k_selection));

  // k is the edit distance in this context

  // generate all combinations of uinque unsigned integers
  prunejuice::pattern::combination::generate_combinations_recursive
    <NCollection, Edge, O, S>
    (n_collection, k_selection, n_items, edge_combinations,
    edge_set_combinations); // sequential routine

  assert(nCk == edge_combinations.size());
  assert(nCk == edge_set_combinations.size());
}

template <typename TemplateGraph, typename TemplateGraphVector, typename S,
  typename TemplateConstraint, typename TemplateConstraintVector>
void generate_k_edit_distance_prototypes_parallel(TemplateGraph& input_template, 
  TemplateGraphVector& prototypes, S& edge_set_combinations, 
  TemplateConstraintVector& template_constraints) {

  {
  #pragma omp parallel for 
  for (size_t i = 0; i < edge_set_combinations.size(); i++) {
    typedef typename TemplateGraph::Vertex_t Vertex_t; 
    typedef typename TemplateGraph::Edge_t Edge_t;
    
    typedef typename TemplateGraph::EdgeListTuple_t EdgeListTuple_t;
    typedef typename TemplateGraph::EdgeListMap_t EdgeListMap_t;
    typedef typename TemplateGraph::EdgeList_t EdgeList_t;    
    typedef typename TemplateGraph::VertexList_t VertexList_t;

    typedef typename TemplateGraph::EdgeListTupleVector_t EdgeListTupleVector_t;
    typedef std::vector<EdgeListTupleVector_t> EdgeListTupleVectors_t; 

    typedef std::unordered_map<Edge_t, std::tuple<Vertex_t, Vertex_t, size_t>> EdgeSet_t; // TODO: ?
    typedef std::vector<EdgeSet_t> EdgeSetVector_t; // TODO: ?

    EdgeListTuple_t new_edge_list = TemplateGraph::filter_edge_list
      //<Edge_t, EdgeListTuple_t, std::unordered_set<Edge_t>> // TODO: define type
      (input_template.edgelist, edge_set_combinations[i]);
    // create a new template graph - a prototype
    TemplateGraph new_template(new_edge_list, input_template.vertex_data, 
      input_template.max_vertex);  
    // Important: must use the max_vertex from input_template.
    // Here the number of vertices is fixed. 

    if (prunejuice::pattern::graph_algorithm::is_connected_component
      <TemplateGraph>(new_template)) {
      // Note: only cycle constraints are generated form the input template,
      // path and TDS constraints are generated for each prototype 
      // (in parallel). template_constraints only contains cycle constraints.

      // identify cycle constraints for the prototype  
      typedef std::unordered_map<size_t, EdgeSet_t> GraphNonLocalPropertiesUnique;
      GraphNonLocalPropertiesUnique prototype_cycles_unique(0);
   
      typedef std::unordered_map<size_t, EdgeSetVector_t> GraphNonLocalPropertiesAllPaths;
      GraphNonLocalPropertiesAllPaths prototype_cycles_all_paths(0);

      GraphNonLocalPropertiesUnique prototype_tds_cycles_unique(0);

      for (size_t j = 0; j < template_constraints.size(); j++) {
        if (template_constraints[j].constraint_type == TemplateConstraint::CYCLE) { 
          if(TemplateGraph::is_subset_edgeset
            (new_template.edgelist_unique, template_constraints[j].edgeset)) {
            assert(j == template_constraints[j].constraint_ID);
            new_template.template_constraints.
              push_back(template_constraints[j].constraint_ID);

            {
            auto find_edgeset_hash = prototype_cycles_unique.
              find(template_constraints[j].edgeset_hash); 
            if (find_edgeset_hash ==  prototype_cycles_unique.end()) {
              prototype_cycles_unique.insert
                ({template_constraints[j].edgeset_hash, 
                template_constraints[j].edgeset});                   
            } else {
              std::cerr << "Error: unexpected item in the map." << std::endl;
              //return; // TODO: graceful exit; return causes OpenMP error   
            }
            }    

            { 
            auto find_edgeset_hash = prototype_cycles_all_paths.
              find(template_constraints[j].edgeset_hash); 
            if (find_edgeset_hash ==  prototype_cycles_all_paths.end()) {
              prototype_cycles_all_paths.insert
                ({template_constraints[j].edgeset_hash, 
                template_constraints[j].edgeset_vector});                   
            } else {
              std::cerr << "Error: unexpected item in the map." << std::endl;
              //return; // TODO: graceful exit; return causes OpenMP error   
            }
            }           
 
          }
        }
      } 

      // TODO: remove template_constraints from TemplateGraph ?     

      // generate prototype constraints 
       
      generate_tds_cycles<Vertex_t, Edge_t>(prototype_cycles_unique, 
        prototype_tds_cycles_unique);
  
      EdgeListTupleVectors_t prototype_tds_cycles_all_paths(0);

      generate_tds_cycle_constraints<Vertex_t, Edge_t, EdgeSet_t,
        GraphNonLocalPropertiesUnique, EdgeListTuple_t, EdgeListTupleVector_t,
        EdgeListTupleVectors_t>(prototype_tds_cycles_unique, 
        prototype_tds_cycles_all_paths);

      EdgeListTuple_t prototype_enumeration_unique = 
        prunejuice::pattern::graphalgorithm::edge_dfs<TemplateGraph>
          (new_template, static_cast<Vertex_t>(0)); 
      
      TemplateConstraintVector protoype_constraints(0);
     
      //populate_template_constraints<Vertex_t, Edge_t, EdgeSet_t, 
      //  TemplateConstraint>(prototype_cycles_unique, prototype_cycles_all_paths, 
      //  prototype_tds_cycles_unique, // TODO: prototype_tds_cycles_all_paths
      //  protoype_constraints);               

      populate_template_constraints<Vertex_t, Edge_t, EdgeSet_t,
        TemplateConstraint, GraphNonLocalPropertiesUnique, 
        GraphNonLocalPropertiesAllPaths, EdgeListTuple_t, 
        EdgeListTupleVector_t, EdgeListTupleVectors_t,
        TemplateConstraintVector>(prototype_cycles_unique, 
        prototype_cycles_all_paths, prototype_tds_cycles_unique, 
        prototype_tds_cycles_all_paths, prototype_enumeration_unique,
        protoype_constraints);

      new_template.template_nonlocal_constraints = std::move(protoype_constraints); // Important:

      // add the template prototype to prototypes 
      {
      #pragma omp critical(prototypes)
      {
        prototypes.push_back(new_template); // TODO: add prototype_template_constraints

      } // #pragma omp citical 
      } 

    } else {
      //std::cout << "Not a connected component." << std::endl;
    }  

  } // for
  } // #pragma omp parallel for
}

template <typename TemplateGraph, typename TemplateGraphVector, 
  typename TemplateConstraint, typename TemplateConstraintVector>
void generate_k_edit_distance_prototypes
  (TemplateGraph& input_template, size_t k_edit_distance, 
  TemplateGraphVector& prototypes, 
  TemplateConstraintVector& template_constraints) {  

  typedef typename TemplateGraph::Edge_t Edge_t;
  typedef std::vector<std::vector<size_t>> O; // indices // TODO
  typedef std::vector<std::unordered_set<Edge_t>> S; // edge_uint // TODO: define this type in template_graph class?

  O edge_combinations(0);
  S edge_set_combinations(0);

  generate_edge_combinations<Edge_t, TemplateGraph, O, S>(input_template, 
    k_edit_distance, edge_combinations, edge_set_combinations);

  generate_k_edit_distance_prototypes_parallel<TemplateGraph, 
    TemplateGraphVector, S, TemplateConstraint, TemplateConstraintVector>
    (input_template, prototypes, edge_set_combinations, template_constraints);
}

template <typename TemplateGraph, typename TemplateGraphVector, 
  typename TemplateGraphVectors, typename TemplateConstraint,
  typename TemplateConstraintVector>
void generate_up_to_k_edit_distance_prototypes 
  (TemplateGraph& input_template, size_t max_k_edit_distance, 
  TemplateGraphVectors& k_edit_distance_prototypes,
  TemplateConstraintVector& template_constraints) {

   assert(k_edit_distance_prototypes.size() > 0); 
   // the first element in k_edit_distance_prototypes is the input template
  
   // TODO: make this loop parallel? // Note: this step can be sequential
   for (size_t k = 1; k <= max_k_edit_distance; k++) {
     generate_k_edit_distance_prototypes<TemplateGraph, TemplateGraphVector, 
       TemplateConstraint, TemplateConstraintVector>
       (input_template, k, k_edit_distance_prototypes[k], template_constraints);      
   } 
}

//------------------------------------------------------------------------------

template <typename Vertex, typename EdgeSet, typename OrderedPathVertices>
void edgeset_to_ordered_path_vertices(EdgeSet& edgeset, 
  OrderedPathVertices& ordered_path_vertices) {
  
  for (auto it = edgeset.begin(); it != edgeset.end(); ++it) {
    //Edge edge_hash = it->first; //ignore
    Vertex s = std::get<0>(it->second);
    Vertex t = std::get<1>(it->second);
    size_t hop = std::get<2>(it->second);
  
    if (hop == 1) {
      ordered_path_vertices[0] = s;
      ordered_path_vertices[hop] = t;
    } else {
      ordered_path_vertices[hop] = t;
    }   
  } 
}



template <typename Vertex, typename EdgeListTuple, typename OrderedPathVertices>
void edgelist_to_ordered_path_vertices(EdgeListTuple& edgelist, 
  OrderedPathVertices& ordered_path_vertices, 
  OrderedPathVertices& ordered_path_indices) {

  std::unordered_map<Vertex, size_t> vertex_path_index_map(0);  
  Vertex current_vertex = 0;
 
  for (size_t i = 0; i < edgelist.size(); i++) {
    Vertex s = std::get<0>(edgelist[i]);
    Vertex t = std::get<1>(edgelist[i]);

    if (i == 0) {
      ordered_path_vertices.push_back(s);

      auto find_s = vertex_path_index_map.find(s);
      if (find_s == vertex_path_index_map.end()) {
        vertex_path_index_map.insert({s, ordered_path_vertices.size() - 1});

        ordered_path_indices.push_back(ordered_path_vertices.size() - 1); 
      }

      ordered_path_vertices.push_back(t);
      current_vertex = t; 

      auto find_t = vertex_path_index_map.find(t);
      if (find_t == vertex_path_index_map.end()) {
        vertex_path_index_map.insert({t, ordered_path_vertices.size() - 1});

        ordered_path_indices.push_back(ordered_path_vertices.size() - 1);
      }

      continue;
    } else {
      if (s != current_vertex) {
        ordered_path_vertices.push_back(s);         
       
        auto find_s = vertex_path_index_map.find(s);
        if (find_s == vertex_path_index_map.end()) {
          vertex_path_index_map.insert({s, ordered_path_vertices.size() - 1});

          ordered_path_indices.push_back(ordered_path_vertices.size() - 1);
        } else {
          ordered_path_indices.push_back(find_s->second);  
        }
           
      }
      ordered_path_vertices.push_back(t);
      current_vertex = t;

      auto find_t = vertex_path_index_map.find(t);
      if (find_t == vertex_path_index_map.end()) {
        vertex_path_index_map.insert({t, ordered_path_vertices.size() - 1});

        ordered_path_indices.push_back(ordered_path_vertices.size() - 1);
      } else {
        ordered_path_indices.push_back(find_t->second); 
      } 

      continue;
    } 
  } // for

  assert(ordered_path_vertices.size() == ordered_path_indices.size()); 
} 

template <typename TemplateGraph, typename TemplateConstraint>
std::string print_cycle_constraint(TemplateConstraint& template_constraint) {
  assert(template_constraint.constraint_type == TemplateConstraint::CYCLE);

  typedef typename TemplateGraph::Vertex_t Vertex_t;  
  typedef std::vector<Vertex_t> OrderedPathVertices;
  std::string constraints = "";
  for (size_t i = 0; i < template_constraint.edgeset_vector.size(); i++) {

    std::string path_string = "";

    size_t path_length = template_constraint.edgeset_vector[i].size() + 1; // Important:   
    OrderedPathVertices ordered_path_vertices(path_length);

    edgeset_to_ordered_path_vertices<Vertex_t>(template_constraint.edgeset_vector[i], 
      ordered_path_vertices);     
    assert(ordered_path_vertices[0] == 
      ordered_path_vertices[ordered_path_vertices.size() - 1]); 

    // path
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {  
      path_string += std::to_string(ordered_path_vertices[j]) + " ";   
    }

    path_string += ": ";

    // path indices 
    // TODO: for now, it assumes unique labels
    // if the path has duplicate intermediate vertices
    // (other than the source), is TDS is 1.
    // Also, need to set the path indices accordingly.  
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      if (j == ordered_path_vertices.size() - 1) {
        path_string += std::to_string(0) + " ";  
      } else {
        path_string += std::to_string(j) + " ";
      }
    }
 
    path_string += ": ";

    // aggregation indices 
    // TODO: for now, it assumes unique labels
    // if the path has duplicate intermediate vertices
    // (other than the source), is TDS is 1.
    // Also, need to set the aggregation vertices accordingly. 
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {    
      path_string += std::to_string(0) + " ";
    }

    // is cyclic : is TDS : do invoke LCC : constraint_ID 
    // TODO: use ConstraintType, set the parameters dynamically 
    // - is cyclic is 1 when source == destination, no need to ack the source
    // - is TDS is 1 if a cycle or a path has duplicate intermediate vertices 
    // (other than the source)
    // TODO: constraint_ID, remove ? 
    // std::string path_parameter_string = ": 1 : 0 : 1 : " + 
    //   std::to_string(template_constraint.constraint_ID);
    std::string path_parameter_string = ": 1 : 0 : 1";  
    //path_string += ": 1 : 0 : 1";
    path_string += path_parameter_string; 

    // write to the file
    constraints += path_string + std::string("\n");
  } // for
  return constraints;
}


template <typename TemplateGraph, typename TemplateConstraint, typename EdgeList>
std::string print_tds_constraint(TemplateConstraint& template_constraint, EdgeList& edge_list) {
  //assert(template_constraint.constraint_type == TemplateConstraint::TDS);
  if (!(template_constraint.constraint_type == TemplateConstraint::TDS || 
    template_constraint.constraint_type == TemplateConstraint::ENUMERATION)) {
    assert(false);
  }

  typedef typename TemplateGraph::Vertex_t Vertex_t;  
  typedef std::vector<Vertex_t> OrderedPathVertices;
  std::string constraints = "";
  for (size_t i = 0; i < template_constraint.edgelist_vector.size(); i++) {

    std::string path_string = "";

    OrderedPathVertices ordered_path_vertices(0);
    OrderedPathVertices ordered_path_indices(0);

    edgelist_to_ordered_path_vertices<Vertex_t>
      (template_constraint.edgelist_vector[i], ordered_path_vertices, 
      ordered_path_indices);
    
    bool is_enum = false;
    if(template_constraint.constraint_type == TemplateConstraint::ENUMERATION){
        // checking whether cover all the edges
        std::vector<std::pair<int, int> > walked_path;
        for (size_t j = 0; j+1 < ordered_path_vertices.size(); j++) {
            if(ordered_path_vertices[j] > ordered_path_vertices[j+1]){
                walked_path.push_back({ordered_path_vertices[j+1], ordered_path_vertices[j]});
            }else{
                walked_path.push_back({ordered_path_vertices[j], ordered_path_vertices[j+1]});
            }
        }
        is_enum = true;
        for(auto e : edge_list){
            bool find = false;
            for(auto s:walked_path){
                if(s.first == e.first && s.second == e.second){
                    find=true;
                    break;
                }
            }
            if(!find){
                is_enum = false;
                break;
            }
        }
        if(!is_enum){
            continue;
        }
    }

    // path
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      path_string += std::to_string(ordered_path_vertices[j]) + " ";
    }

    path_string += ": ";

    // path indices
    for (size_t j = 0; j < ordered_path_indices.size(); j++) {
      path_string += std::to_string(ordered_path_indices[j]) + " ";
    }

    path_string += ": ";

    // aggregation indices // TODO: set the aggregation vertices 
    for (size_t j = 0; j < ordered_path_vertices.size(); j++) {
      path_string += std::to_string(1) + " ";
    }

    path_string += ": ";

    bool is_cyclic = ordered_path_vertices[0] == 
      ordered_path_vertices[ordered_path_vertices.size() - 1] ? true : false;   

    path_string += std::to_string(is_cyclic) + " "; 

    path_string += ": ";

    bool is_TDS = true; // TODO: constraint type ?

    path_string += std::to_string(is_TDS) + " ";

    path_string += ": "; 

    bool do_lcc = true; // TODO: ?

    path_string += std::to_string(do_lcc) + " "; 

    // path_string += ": ";

    // path_string +=  std::to_string(template_constraint.constraint_ID); // TODO: constraint_ID, remove ?
    if(template_constraint.constraint_type == TemplateConstraint::ENUMERATION){
        // constraints += path_string+std::string("\n");
    //   if(i == template_constraint.edgelist_vector.size()-1){
        return path_string+std::string("\n");
    //   }
    }else{
        constraints += path_string+std::string("\n");
    }
  } // for
    return constraints;
}


/**
 * DFS walk on the TDS constraint
 */
template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList, typename OrderedPath, typename Visited>
void dfs_recursive(VertexList& vertices, Vertex vertex_count,
  VertexList& vertex_degree, EdgeList& edges,
  Vertex source_vertex, Vertex v, OrderedPath& walk_history, 
  Visited& visited, size_t r) {

  visited[static_cast<size_t>(v)] = 1;  

  for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
    Vertex v_nbr = edges[e];

    if (visited[static_cast<size_t>(v_nbr)] == 0) { 
      dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
        (vertices, vertex_count, vertex_degree, edges, source_vertex, v_nbr, 
        walk_history, visited, r + 1);

    } else {
      // std::cout << "[" << v << " - " << v_nbr << "], "; // Test
    }
    
  } // for 
}

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList>
void dfs_parallel(VertexList& vertices, Vertex vertex_count,
  VertexList& vertex_degree, EdgeList& edges) {
  
  {
  //#pragma omp parallel for schedule(static, 1)
  for (Vertex v = 0; v < vertex_count; v++) 
  {

    std::vector<uint8_t> visited(vertex_count);
    size_t r = 1; // walk step, initialized to 1

    typedef std::unordered_map<Vertex, size_t> OrderedPath;
    OrderedPath walk_history(0);      

    auto find_v = walk_history.find(v); 
    if (find_v == walk_history.end()) {
      walk_history.insert({v, r});  
    } else {
      std::cerr << "Error: unexpected item in the map." << std::endl;
    }

    dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
      (vertices, vertex_count, vertex_degree, edges, v, v, walk_history, 
      visited, r);
  } // for
  } 
} 
 
template <typename Vertex, typename Edge, typename VertexList,
  typename EdgeList>
void dfs_single_source(VertexList& vertices, Vertex vertex_count,
  VertexList& vertex_degree, EdgeList& edges, Vertex v) {
  
  {
  //#pragma omp parallel for schedule(static, 1)
  for (Vertex v = 0; v < vertex_count; v++) 
  {

    std::vector<uint8_t> visited(vertex_count);
    if (vertex_degree[v] < 1) {
      continue;
    }

    size_t r = 1; // walk step, initialized to 1

    typedef std::unordered_map<Vertex, size_t> OrderedPath;
    OrderedPath walk_history(0);      

    auto find_v = walk_history.find(v); 
    if (find_v == walk_history.end()) {
      walk_history.insert({v, r});  
    } else {
      std::cerr << "Error: unexpected item in the map." << std::endl;
    }

    dfs_recursive<Vertex, Edge, VertexList, EdgeList, OrderedPath>
      (vertices, vertex_count, vertex_degree, edges, v, v, walk_history, 
      visited, r);
    break; // only require one successful walk from one source
  } // for
  } 
} 

////////////////////////////////////////////////////////////////////////////////

/**
 * This routine is called by multiple graphs is parallel, so the routine itself
 * is siquential 
 */
template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList> 
bool is_connected_component(VertexList& vertices, Vertex vertex_count, 
  EdgeList& edges) {

  std::vector<uint8_t> visited(vertex_count);
  for (auto& v : visited) {
    v = 0;
  }

  Vertex source = 0;
  visited[static_cast<size_t>(source)] = 1;  
  bool finished = false;
  bool is_cc = true;

  // TODO: just count the number of visited vertices?

  // BFS 
  do {
    finished = true;
    for (Vertex v = 0; v < vertex_count; v++) {
      if (visited[static_cast<size_t>(v)] == 1) { 
        for (Edge e = vertices[v]; e < vertices[v + 1]; e++) {
          Vertex v_nbr = edges[e];  
          if (visited[static_cast<size_t>(v_nbr)] == 0) {  
            visited[static_cast<size_t>(v_nbr)] = 1;
            finished = false; 
          } // if
        } // for 
      } // if  
    } // for
  } while(!finished);

  for (Vertex v = 0; v < vertex_count; v++) {        
    if (visited[static_cast<size_t>(v)] != 1) {
      is_cc = false; // the input graph is not a connected component 
      break;
    }  
  }

  return is_cc;
}  

////////////////////////////////////////////////////////////////////////////////

/**
 * Generate combinations
 */
//template <typename N, typename K, typename NCK>
template <typename NCollection, typename ItemType, typename O, 
  typename CombinationBuffer, typename S>
void generate_combinations_recursive(size_t s, size_t n, size_t k, size_t r, 
  CombinationBuffer& buffer, O& combinations, NCollection& n_items, 
  S& combinations_items) {
  if (r == k) {
    combinations.push_back(buffer);
    
    std::unordered_set<ItemType> buffer_set(0); 
    for (auto i : buffer) {
      buffer_set.insert(n_items[i]);    
    }
    assert(buffer.size() == buffer_set.size());
    combinations_items.push_back(buffer_set);

    return; 
  }

  for (size_t i = s; i < n ; i++) {
    buffer[r] = i;
    generate_combinations_recursive<NCollection, ItemType, O, CombinationBuffer, S>
      (i + 1, n, k, r + 1, buffer, combinations, n_items, combinations_items);
  } // for  
}

//template <typename N, typename K, typename NCK>
template <typename NCollection, typename ItemType, typename O, typename S>
void generate_combinations_recursive(size_t n, size_t k, 
  NCollection& n_items, O& combinations, S& combinations_items) {

  typedef std::vector<size_t> CombinationBuffer;
  CombinationBuffer buffer(k); // size_t - array index type

  generate_combinations_recursive<NCollection, ItemType, O, CombinationBuffer, S>
    (0, n, k, 0, buffer, combinations, n_items, combinations_items);
    // size_t s, size_t n, size_t k, size_t r, ...	
}

template <typename N, typename K, typename NCK, 
  typename CombinationBitSet, typename CombinationCollection>
void generate_combintions_parallel_3(N _n, K _k, NCK _nCk, 
  CombinationCollection & combinations) { 
  size_t n = static_cast<size_t>(_n);
  size_t k = static_cast<size_t>(_k); 
  size_t nCk = static_cast<size_t>(_nCk);

  assert(nCk == combinations.size()); 
} 

template <typename N, typename K, typename NCK, 
  typename CombinationBitSet, typename CombinationCollection>
void generate_combintions_parallel_2(N _n, K _k, NCK _nCk, 
  CombinationCollection & combinations) {
  size_t n = static_cast<size_t>(_n);
  size_t k = static_cast<size_t>(_k); 
  size_t nCk = static_cast<size_t>(_nCk);

  size_t s = n - k; // number of items to remove 

  assert(nCk == combinations.size()); 

  {
  //#pragma omp parallel for schedule(static, 1)
  for (size_t i = 0; i < nCk; ++i) {

    size_t thread_ID = omp_get_thread_num();
    
    for (size_t j = i*s; j < i*s + s; j++) {
      size_t r = j % n;
      
      assert(r < n);    
      assert(combinations[i].test(r) == false);
      combinations[i].set(r);
    }
  }
 
  } 
   
} 

template <typename N, typename K, typename NCK, 
  typename CombinationBitSet, typename CombinationCollection>
void generate_combintions_parallel(N _n, K _k, NCK _nCk, 
  CombinationCollection & combinations) {
  size_t n = static_cast<size_t>(_n);
  size_t k = static_cast<size_t>(_k); 
  size_t nCk = static_cast<size_t>(_nCk);

  assert(nCk == combinations.size()); 

  size_t s = n - k; // number of items to remove 

  {
  #pragma omp parallel for schedule(static, 1)
  for (size_t i = 0; i < nCk; ++i) {

    size_t thread_ID = omp_get_thread_num();
    
    size_t r = i % n; // r is always smaller than n
    assert(r < n);

    size_t offset = i / n; // integer division 
    assert(offset < n);
    // r is the first item
    assert(combinations[i].test(r) == false); 
    combinations[i].set(r);
    
    for (size_t j = 0; j < s - 1; j++) {

      size_t offset_j = offset + j + 1;

      if (offset_j > (n - 1)) {
        offset_j = offset_j - n;
        assert(offset_j < n);
      }  

      size_t next_item = r + offset_j; 
      // next_item can be > n
 
      if (next_item > (n - 1)) {
        next_item = next_item - n;
        assert (next_item < n);  
      }
      combinations[i].set(next_item);    
    } // for

 
  } // for 
 
  } 
   
} 

////////////////////////////////////////////////////////////////////////////////

/**
 * Graph construction
 */
template <typename EdgeListTuple>
void write_edge_list_file(EdgeListTuple& edge_list, size_t unique_ID, 
  std::string output_filepath) {
  std::string output_filename = output_filepath + "/edgelist_" + 
    std::to_string(unique_ID);
  std::ofstream output_file(output_filename, std::ofstream::out);
  for (size_t e = 0; e < edge_list.size(); e++) {    
    output_file << std::get<0>(edge_list[e]) << " " 
      << std::get<1>(edge_list[e]) << "\n";  
  }  
  output_file.close();
} 

template <typename Vertex, typename VertexData, typename VertexDataList>
void read_vertex_data_file(const std::string vertex_data_input_filename, 
  VertexDataList& vertex_data) {
  std::ifstream vertex_data_input_file(vertex_data_input_filename,
    std::ifstream::in);
  std::string line;
  while (std::getline(vertex_data_input_file, line)) {
    std::istringstream iss(line);
    Vertex v_source(0);
    VertexData v_data(0);
    iss >> v_source >> v_data;
    vertex_data.push_back(v_data);
  }
  vertex_data_input_file.close();
}

template <typename Vertex, typename Edge, typename EdgeListTuple>
Edge read_edge_list_file(const std::string input_filename, EdgeListTuple& edge_list,
  Vertex& max_vertex, Edge skip_lines_count = 0) {
  std::ifstream input_file(input_filename, std::ifstream::in);
  Edge edge_count(0);
  Vertex s(0), t(0);
  std::string line;

  while(std::getline(input_file, line)) {
    std::istringstream iss(line);
    edge_count++;
    //if (edge_count >= skip_lines_count || skip_lines_count == 0) {
      iss >> s >> t;
      edge_list.push_back(std::forward_as_tuple(s, t));

      auto tmp_max_vertex = s >= t ? s : t;

      if (max_vertex < tmp_max_vertex) {
        max_vertex = tmp_max_vertex;
      }

    //}
  }
  input_file.close();
  edge_count-= skip_lines_count;
  assert(edge_count > 0);
  return edge_count;
}

template <typename EdgeListTuple, typename EdgeList>
void generate_edge_list(EdgeListTuple& edge_list, EdgeList& edges) {
  for (size_t e = 0; e < edge_list.size(); e++) {
    edges.push_back(std::get<1>(edge_list[e]));
  }
}

template <typename EdgeListTuple>
EdgeListTuple directed_to_undirected_edge_list
  (EdgeListTuple& directed_edge_list) {  
  EdgeListTuple undirected_edge_list(directed_edge_list);
  assert(directed_edge_list.size() == undirected_edge_list.size());
  for (size_t e = 0; e < directed_edge_list.size(); e++) {
    undirected_edge_list.
      push_back(std::forward_as_tuple(std::get<1>(directed_edge_list[e]), 
      std::get<0>(directed_edge_list[e]))); 
  } // for
  return undirected_edge_list; 
} 

template <typename Vertex, typename Edge, typename EdgeListTuple, 
  typename EdgeListMap>
void generate_unique_edge_list(EdgeListTuple& edge_list, 
  EdgeListMap& edge_list_unique) {

  for (auto edge : edge_list) {
    auto s = std::get<0>(edge);
    auto t = std::get<1>(edge);
    // unique edge identifier
    std::bitset<VERTEX_COUNT> edge_bitset; // TODO: this is actually vertex bitset / edge hash? 
    edge_bitset.set(static_cast<size_t>(s));
    edge_bitset.set(static_cast<size_t>(t));      
 
    Edge edge_uint = static_cast<Edge>(edge_bitset.to_ullong());

    edge_list_unique.insert({edge_uint, edge});
  }      
}  

// retun a new edge list without the edges in edge_list_filter 
template <typename Edge, typename EdgeListTuple, typename EdgeListFilter>
EdgeListTuple filter_edge_list(EdgeListTuple& edge_list, 
  EdgeListFilter& edge_list_filter) {

  EdgeListTuple new_edge_list(0);

  for (auto edge : edge_list) {
    auto s = std::get<0>(edge);
    auto t = std::get<1>(edge);
   
    // unique edge identifier
    std::bitset<VERTEX_COUNT> edge_bitset;
    edge_bitset.set(static_cast<size_t>(s));
    edge_bitset.set(static_cast<size_t>(t));

    Edge edge_uint = static_cast<Edge>(edge_bitset.to_ullong());

    auto find_edge = edge_list_filter.find(edge_uint); // ignore if found
    if (find_edge == edge_list_filter.end()) { 
      // edge not found, add it to the new_edge_list
      new_edge_list.push_back(edge);   
    }   
  }
  
  return new_edge_list; 
}   

template <typename EdgeListTuple, typename VertexList, typename Vertex>          
Vertex generate_vertex_list(EdgeListTuple& edge_list, VertexList& vertices,
  VertexList& vertex_degree, const Vertex max_vertex) {

  //std::ofstream vertex_file("vertex_file_tmp", std::ofstream::out);

  Vertex vertex_count = 0;
  //Vertex max_vertex = std::get<0>(edge_list[edge_list.size() - 1]);
  Vertex l = 0; // edge list index
  Vertex degree = 0;
  Vertex current_vertex = vertex_count;
  Vertex source;
  //Vertex target;

  do {
    auto edge = edge_list[l];
    source = std::get<0>(edge);
    //target = std::get<1>(edge);
    if (source == current_vertex) {
      degree++;
      l++;
    } else {
      vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                     vertices[vertices.size() - 1]));
      vertex_degree.push_back(degree);
      // update vertices array
      degree = 0;
      vertex_count++;
      current_vertex = vertex_count;
    }
  } while(current_vertex <= max_vertex);

  // add the last dummy vertex to the vertex list
  vertices.push_back(vertices.size() == 0 ? 0 : (vertex_degree[vertex_degree.size() - 1] +
                                                 vertices[vertices.size() - 1]));

  return vertex_count;
}

/**
 * Generate CSR graph from a given edgelist
 */
// TODO: this function assumes the input edgelist is undirected? 
// Note: it also works for a directed edgelist
template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeListTuple, typename EdgeListMap>
void generate_graph(Vertex& vertex_count, VertexList& vertices, 
  VertexList& vertex_degree, EdgeList& edges, EdgeListTuple& edge_list,
  EdgeListMap& edge_list_unique, Vertex& max_vertex) {
  
  std::stable_sort(edge_list.begin(), edge_list.end(),
    [](const std::tuple<Vertex, Vertex>& a,
       const std::tuple<Vertex, Vertex>& b) -> bool {
         return std::get<0>(a) < std::get<0>(b);
       });
  // TODO: make sure this is not the issue again 
  // Important: max_vertex > 0 means, max_vertex is defined already
  if (max_vertex < 1) {  
    for (auto& e : edge_list) {
      auto s = std::get<0>(e);
      auto t = std::get<1>(e);
  
      auto tmp_max_vertex = s >= t ? s : t;

      if (max_vertex < tmp_max_vertex) {
          max_vertex = tmp_max_vertex;
      }
    } 
  } // if

  vertex_count = generate_vertex_list(edge_list, vertices, vertex_degree, 
    max_vertex);

  // {
  // #pragma omp parallel for
  for (size_t v = 0; v < vertices.size() - 1; v++) {
    size_t start = vertices[v];
    size_t end = vertices[v + 1];
    std::stable_sort(edge_list.begin() + start,
                     edge_list.begin() + end,
       [](const std::tuple<Vertex, Vertex>& a,
          const std::tuple<Vertex, Vertex>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });
  }
  // }

  // generate edge list
  generate_edge_list(edge_list, edges);

  if (edge_list_unique.size() < 1) {
    generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
      (edge_list, edge_list_unique);
  }        
}

template <typename Vertex, typename Edge, typename VertexList, 
  typename EdgeList, typename EdgeListTuple>
void generate_graph(Vertex& vertex_count, VertexList& vertices, 
  VertexList& vertex_degree, EdgeList& edges, EdgeListTuple& edge_list,
  Vertex& max_vertex) {

  // TODO: imporve
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> T;  
  T dummy_edge_list_unique;
  dummy_edge_list_unique.insert({0, std::forward_as_tuple(0, 0)});
 
  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple, T>
    (vertex_count, vertices, vertex_degree, edges, edge_list, 
    dummy_edge_list_unique, max_vertex); 
} 

////////////////////////////////////////////////////////////////////////////////
template <typename Vertex, typename Edge, typename EdgeListTuple, typename VertexData, typename VertexDataList>
Edge read_query(std::string input_filename, EdgeListTuple& edge_list,
  Vertex& max_vertex, VertexDataList& vertex_data) {
  std::ifstream input_file(input_filename, std::ifstream::in);
  Edge edge_count(0);
  Vertex s(0), t(0);
  std::string line;
  while(std::getline(input_file, line)) {
    std::istringstream iss(line);
    std::string none;
    
    if(line[0] == 'e'){
      edge_count+=2;
      iss>>none>>s>>t;
      edge_list.push_back(std::forward_as_tuple(s, t));
      edge_list.push_back(std::forward_as_tuple(t, s));
      auto tmp_max_vertex = s >= t ? s : t;
      if (max_vertex < tmp_max_vertex) {
        max_vertex = tmp_max_vertex;
      }
    }else if(line[0] == 'v'){
      Vertex v_source(0);
      VertexData v_data(0);
      iss>>none>>v_source>>v_data;
      vertex_data.push_back(v_data);
    }
  }
  input_file.close();
  return edge_count;
}


std::string generate_constraints(std::string filename){
  typedef uint64_t Vertex;
  typedef uint64_t Edge;
  typedef uint64_t VertexData;
  typedef uint64_t EdgeData;

  typedef std::vector<std::tuple<Vertex, Vertex>> EdgeListTuple;
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> EdgeListMap;
  typedef std::vector<Vertex> EdgeList;
  typedef std::vector<Edge> VertexList;
  typedef std::vector<VertexData> VertexDataList;

  typedef std::vector<EdgeListTuple> EdgeListTupleVector;
  typedef std::vector<EdgeListTupleVector> EdgeListTupleVectors;

  typedef prunejuice::pattern::template_graph<Vertex, Edge, VertexData, 
    VertexList, EdgeList, EdgeListTuple, EdgeListMap, VertexDataList> TemplateGraph;
  typedef std::vector<TemplateGraph> TemplateGraphVector;
  typedef std::vector<TemplateGraphVector> TemplateGraphVectors; 

  // edge hash, s, t, hop ID - edge hash is used by cycles, hop ID is used by cycles
  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex, size_t>> EdgeSet;
  typedef std::vector<EdgeSet> EdgeSetVector;

  typedef prunejuice::pattern::template_constraint<Vertex, Edge, EdgeSet, 
    EdgeSetVector, EdgeListTupleVector> TemplateConstraint;
  typedef std::vector<TemplateConstraint> TemplateConstraintVector;

  typedef prunejuice::pattern::file_utilities<TemplateGraph> FileUtilities;

  typedef std::unordered_map<VertexData, std::unordered_map<Edge, 
    std::tuple<Vertex, Vertex>>> VetexDataVertexPairs;

  //////////////////////////////////////////////////////////////////////////////

  // graph construction
  EdgeListTuple edge_list(0);
  EdgeListMap edge_list_unique(0);
  EdgeList edges(0);
  VertexList vertices(0);
  VertexList vertex_degree(0);
  VertexDataList vertex_data(0);

  Edge edge_count = 0;
  Vertex vertex_count = 0;
  Vertex max_vertex = 0;
  Edge max_degree = 0;

  std::chrono::time_point<std::chrono::steady_clock> global_start_time =
    std::chrono::steady_clock::now();

  edge_count = read_query<Vertex, Edge, EdgeListTuple, VertexData, VertexDataList>(filename, edge_list,
    max_vertex, vertex_data);

  std::stable_sort(edge_list.begin(), edge_list.end(),
    [](const std::tuple<Vertex, Vertex>& a,
       const std::tuple<Vertex, Vertex>& b) -> bool {
         return std::get<0>(a) < std::get<0>(b);
       });

  // generate vetex list
  vertex_count = generate_vertex_list(edge_list, vertices, vertex_degree, max_vertex);

  for (size_t v = 0; v < vertices.size() - 1; v++) {
    size_t start = vertices[v];
    size_t end = vertices[v + 1];
    std::stable_sort(edge_list.begin() + start,
                     edge_list.begin() + end,
       [](const std::tuple<Vertex, Vertex>& a,
          const std::tuple<Vertex, Vertex>& b) -> bool {
            return std::get<1>(a) < std::get<1>(b);
          });
  }

  // generate edge list
  generate_edge_list(edge_list, edges);
  generate_unique_edge_list<Vertex, Edge, EdgeListTuple, EdgeListMap>
   (edge_list, edge_list_unique);

  // read vertex data
  TemplateGraph input_template(edge_list, vertex_data, 0);  

  //////////////////////////////////////////////////////////////////////////////
  
  // build in-memory CSR graph (directed) // Note: not used
  global_start_time =
    std::chrono::steady_clock::now();

  EdgeListTuple edge_list_directed(0);
  EdgeList edges_directed(0);
  VertexList vertices_directed(0);
  VertexList vertex_degree_directed(0);

  Vertex vertex_count_directed = 0;
  Vertex max_vertex_directed = 0; 

  for (auto& e : edge_list_unique) {
    edge_list_directed.push_back(std::forward_as_tuple(std::get<0>(e.second), 
      std::get<1>(e.second)));
  }

  generate_graph<Vertex, Edge, VertexList, EdgeList, EdgeListTuple>
    (vertex_count_directed, vertices_directed, vertex_degree_directed, 
    edges_directed, edge_list_directed, max_vertex_directed);

  typedef std::vector< std::unordered_map<size_t, std::vector<std::tuple<Vertex, Vertex>> > > 
    VertexNonLocalProperties;

  // generating cycles
  VertexNonLocalProperties vertex_cycles(vertex_count);

  typedef std::vector< std::unordered_map<size_t, EdgeSet> > VertexNonLocalPropertiesUnique;
  VertexNonLocalPropertiesUnique vertex_cycles_unique(vertex_count);
  find_cycles_parallel<Vertex, Edge, VertexList, EdgeList>(vertices, 
    vertex_count, vertex_degree, edges, vertex_cycles, vertex_cycles_unique);

  // vector based path
  typedef std::unordered_map<size_t, std::vector<std::tuple<Vertex, Vertex>> > 
    NonLocalProperties; 
  NonLocalProperties graph_cycles(0);

  // identify unique cycles // TODO: graph_cycles is not used, remove? 
  for (size_t i = 0; i < vertex_cycles.size(); i++) {
    for (auto& j : vertex_cycles[i]) { 
      auto find_path = graph_cycles.find(j.first);
      if (find_path == graph_cycles.end()) {
        graph_cycles.insert({j.first, j.second});   
      }    
    } // for       
  } // for 
  
  typedef std::unordered_map<size_t, EdgeSet> GraphNonLocalPropertiesUnique; 
  GraphNonLocalPropertiesUnique graph_cycles_unique(0);

  typedef std::unordered_map<size_t, EdgeSetVector> GraphNonLocalPropertiesAllPaths;
  GraphNonLocalPropertiesAllPaths graph_cycles_all_paths(0); // grouped by identical path / constraint 

  find_unique_cycles<Vertex, Edge>(vertex_cycles_unique, graph_cycles_unique, 
    graph_cycles_all_paths);


  // TODO: move it to a function 
  for (auto i : input_template.vertex_data_vertex_pairs) {
    for (auto j : i.second) {
      Vertex s = std::get<0>(j.second);
      Vertex t = std::get<1>(j.second);   

      EdgeSet forward_path = TemplateGraph::edgelisttuple_to_edgeset
        ( prunejuice::pattern::graph_algorithm::get_shortest_path
        <Vertex, Edge, EdgeListTuple, TemplateGraph>(input_template, s, t) );  

      EdgeSet reverse_path = TemplateGraph::edgelisttuple_to_edgeset
        ( prunejuice::pattern::graph_algorithm::get_shortest_path
        <Vertex, Edge, EdgeListTuple, TemplateGraph>(input_template, t, s) );

    } // for 
  } // for 

  typedef std::unordered_map<Edge, std::tuple<Vertex, Vertex>> TDSEdgeSet;  
  TDSEdgeSet tds_cycles_edge_set;  
  generate_tds_cycles<Vertex, Edge>(graph_cycles_unique, tds_cycles_edge_set); // single edgelist 
  // Note: not used


  GraphNonLocalPropertiesUnique graph_tds_cycles_unique(0); // map of edge sets
 
  generate_tds_cycles<Vertex, Edge>(graph_cycles_unique, graph_tds_cycles_unique);

  EdgeListTupleVectors graph_tds_cycles_all_paths(0);  

  generate_tds_cycle_constraints<Vertex, Edge, EdgeSet, 
    GraphNonLocalPropertiesUnique, EdgeListTuple, EdgeListTupleVector, 
    EdgeListTupleVectors>(graph_tds_cycles_unique, graph_tds_cycles_all_paths, 
    true);

  // generate enumeration constraint
  EdgeListTuple graph_enumeration_unique = 
    prunejuice::pattern::graphalgorithm::edge_dfs<TemplateGraph>
      (input_template, static_cast<Vertex>(0));
  
  // Test
  for (Vertex v = 0; v < input_template.vertex_count; v++) {
    EdgeListTuple graph_enumeration_all_paths = 
      prunejuice::pattern::graphalgorithm::edge_dfs<TemplateGraph>
        (input_template, v);
  }

  TemplateConstraintVector template_constraints(0);  
  TemplateConstraintVector cycle_constraints(0);

  populate_template_constraints<Vertex, Edge, EdgeSet, TemplateConstraint, 
    GraphNonLocalPropertiesUnique, GraphNonLocalPropertiesAllPaths, 
    EdgeListTuple, EdgeListTupleVector, EdgeListTupleVectors,
    TemplateConstraintVector>
    (graph_cycles_unique, graph_cycles_all_paths, graph_tds_cycles_unique, 
    graph_tds_cycles_all_paths, graph_enumeration_unique, 
    template_constraints); 

  // Note: only cycle constraints are generated form the input template,
  // path and TDS constraints are generated for each prototype (in parallel).
  populate_cycle_constraints<Vertex, Edge, EdgeSet, TemplateConstraint>
    (graph_cycles_unique, graph_cycles_all_paths, cycle_constraints);

  //---------------------------------------------------------------------------- 

  // input template constraints  

  input_template.template_nonlocal_constraints = template_constraints; // Important: 

  // identify template constraints for the input_template 
  for (size_t j = 0; j < template_constraints.size(); j++) {
    if (TemplateGraph::is_subset_edgeset
      (input_template.edgelist_unique, template_constraints[j].edgeset)) {
      assert(j == template_constraints[j].constraint_ID);
      input_template.template_constraints.
        push_back(template_constraints[j].constraint_ID);
    }
  } // for  

  //////////////////////////////////////////////////////////////////////////////
  // generate up to k-edit distance template prototypes
  size_t k_input = 1;
  TemplateGraphVectors k_edit_distance_prototypes(k_input + 1); 

  // the first element in k_edit_distance_prototypes is the input template
  k_edit_distance_prototypes[0].push_back(input_template);
//#ifdef ENABLE_BLOCK
  generate_up_to_k_edit_distance_prototypes<TemplateGraph, TemplateGraphVector, 
    TemplateGraphVectors, TemplateConstraint, TemplateConstraintVector>
    (input_template, k_input, k_edit_distance_prototypes, cycle_constraints); //template_constraints);
//#endif
  // Test
  //////////////////////////////////////////////////////////////////////////////
 
  // create input and output directories / write template prototypes to files
  // create_template_prototype_directories<TemplateGraphVectors, FileUtilities>
  //   (k_edit_distance_prototypes, prototype_dir_name);
    std::vector<std::pair<int, int> > edge_pair;
    for(auto v : edge_list_directed){
        edge_pair.push_back({std::get<0>(v), std::get<1>(v)});
    }

  // outputing the constraints
  std::string constraints = "";
  for (size_t k = 0;  k < k_edit_distance_prototypes.size(); k++) {
    if (k_edit_distance_prototypes[k].size() < 1) {
      continue;
    }
    // k = k_edit_distance_prototypes.size()-1;
    for(size_t p=0; p<k_edit_distance_prototypes[k].size(); p++){
        auto template_graph = k_edit_distance_prototypes[k][p]; // pick the 1-st template_graph in default
        for (size_t i = 0; i < template_graph.template_nonlocal_constraints.size(); i++) {
#ifdef ENABLE_CONSTRAINT
            if (template_graph.template_nonlocal_constraints[i].constraint_type == TemplateConstraint::CYCLE) {
                constraints += print_cycle_constraint<TemplateGraph, TemplateConstraint>(template_graph.template_nonlocal_constraints[i]);
            }
#endif
#ifdef ENABLE_ENUMERATION
            if (template_graph.template_nonlocal_constraints[i].constraint_type == TemplateConstraint::ENUMERATION) {
                std::string enum_tds = print_tds_constraint<TemplateGraph, TemplateConstraint>(template_graph.template_nonlocal_constraints[i], edge_pair);
                constraints += enum_tds;
                if(enum_tds.size() != 0){
                    goto out;
                }
            }
#endif
        }
    }
    // break;
  }
  out:
  std::string result(constraints, 0, constraints.size()-1);
  return result;
}
