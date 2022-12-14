// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_DELEGATE_PARTITIONED_GRAPH_HPP_INCLUDED
#define HAVOQGT_MPI_DELEGATE_PARTITIONED_GRAPH_HPP_INCLUDED

#include <limits>
#include <utility>
#include <stdint.h>
#include <functional>

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/map.hpp>
#include <boost/interprocess/containers/map.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

#include <havoqgt/mpi.hpp>
#include <havoqgt/utilities.hpp>
// #include <havoqgt/cache_utilities.hpp>
#include <havoqgt/detail/iterator.hpp>
#include <havoqgt/impl/edge_partitioner.hpp>
#include <havoqgt/impl/edge_node_identifier.hpp>

#ifdef DEBUG_DPG
 #warning Debug MACRO is for delegate_partitioned_graph.
 #define IS_DEBUGING true
#else
 #define IS_DEBUGING false
#endif

namespace havoqgt {

namespace bip = boost::interprocess;

/**
 * Delegate partitioned graph using MPI for communication.
 *
 * @todo Test using simple deterministic patterns.
 * @todo Add edge_data
 * @todo Make vertex_iterator a random access iterator
 * @todo Add invalid bit or state to vertex_locator
 * @todo Verify low-degree CSR creation:  ipp line 167
 * @todo Boostify controller locator
 */
template <typename Allocator>
class delegate_partitioned_graph {
 public:


  //////////////////////////////////////////////////////////////////////////////
  // Class Objects
  //////////////////////////////////////////////////////////////////////////////
  /// Object that uniquely locates a vertex, both MPI rank and local offset.
  class vertex_locator;
  /// Edge Iterator class for delegate partitioned graph
  class edge_iterator;
  /// Vertex Iterator class for delegate partitioned graph
  class vertex_iterator;
  /// Vertex Data storage
  template <typename T, typename AllocatorOther>
  class vertex_data;
  /// Edge Data storage
  template <typename T, typename AllocatorOther>
  class edge_data;
  /// Stores information about owned vertices
  class vert_info;

  enum ConstructionState { New, MetaDataGenerated, EdgeStorageAllocated,
    LowEdgesPartitioned, HighEdgesPartitioned, GraphReady};

  //////////////////////////////////////////////////////////////////////////////
  // Public Member Functions
  //////////////////////////////////////////////////////////////////////////////
  /// Constructor that initializes given and unsorted sequence of edges

  template <typename Container>
  delegate_partitioned_graph(const Allocator& allocator,
                             MPI_Comm mpi_comm,
                             Container& edges, 
                             uint64_t max_vertex,
                             uint64_t delegate_degree_threshold,
                             uint64_t _node_partitions,
                             uint64_t chunk_size,
                             ConstructionState stop_after = GraphReady
                             );

  template <typename Container, typename edge_data_type>  
  delegate_partitioned_graph(const Allocator& allocator,
                             MPI_Comm mpi_comm,
                             Container& edges, 
                             uint64_t max_vertex,
                             uint64_t delegate_degree_threshold,
                             uint64_t _node_partitions,
                             uint64_t chunk_size,
                             edge_data_type& _edge_data,  
                             bool _has_edge_data = true,
                             ConstructionState stop_after = GraphReady
                             );

  template <typename Container, typename edge_data_type>
  void complete_construction(const Allocator& allocator,
    MPI_Comm mpi_comm, Container& edges, edge_data_type& _edge_data);
  void print_graph_statistics();

  /// Converts a vertex_locator to the vertex label
  uint64_t locator_to_label(vertex_locator locator) const;

  /// Converts a vertex label to a vertex_locator
  vertex_locator label_to_locator(uint64_t label) const;

  /// Returns a begin iterator for edges of a vertex
  edge_iterator edges_begin(vertex_locator locator) const;

  /// Returns an end iterator for edges of a vertex
  edge_iterator edges_end(vertex_locator locator) const;

  /// Returns the degree of a vertex
  uint64_t degree(vertex_locator locator) const;

  /// Returns the local degree of a vertex
  uint64_t local_degree(vertex_locator locator) const;

  /// Returns the outgoing degree of a vertex
  uint64_t outgoing_degree(vertex_locator locator) const;

  /// Returns the incoming degree of a vertex
  u_int64_t incoming_degree(vertex_locator locator) const;

  /// Returns a begin iterator for all local vertices
  vertex_iterator vertices_begin() const;

  /// Returns an end iterator for all local vertices
  vertex_iterator vertices_end() const;

  /// Returns a begin iterator for all delegate vertices
  vertex_iterator delegate_vertices_begin() const;

  /// Returns an end iterator for all delegate vertices
  vertex_iterator delegate_vertices_end() const;

  /// Tests if vertex label is a delegate
  bool is_label_delegate(uint64_t label) const;

  /// Creates vertex_data of type T
  template <typename T, typename AllocatorOther>
  vertex_data<T, AllocatorOther>* create_vertex_data(
      const AllocatorOther&,
      const char *obj_name = nullptr) const;

  /// Creates vertex_data of type T, with initial value
  template <typename T, typename AllocatorOther>
  vertex_data<T, AllocatorOther>* create_vertex_data(
      const T& init, const AllocatorOther&,
      const char *obj_name = nullptr) const;

  /// Creates edge_data of type T
  template <typename T, typename AllocatorOther>
  edge_data<T, AllocatorOther>* create_edge_data(
      const AllocatorOther&,
      const char *obj_name = nullptr) const;

  /// Creates edge_data of type T, with initial value
  template <typename T, typename AllocatorOther>
  edge_data<T, AllocatorOther>* create_edge_data(
      const T& init, const AllocatorOther&,
      const char *obj_name = nullptr) const;

  size_t num_local_vertices() const {
    return m_owned_info.size();
  }

  uint64_t max_global_vertex_id() {
    return m_global_max_vertex;
  }

  uint64_t max_local_vertex_id() {
    return m_max_vertex;
  }

  size_t num_delegates() const {
    return m_delegate_degree.size();
  }

  uint32_t master(const vertex_locator& locator) const {
    return locator.m_local_id % m_mpi_size;
  }

  typedef typename bip::vector<vertex_locator, other_allocator<Allocator, vertex_locator> >
      ::const_iterator controller_iterator;

  controller_iterator controller_begin() const {
    return m_controller_locators.begin();
  }

  controller_iterator controller_end()   const {
    return m_controller_locators.end();
  }

  bool compare(delegate_partitioned_graph<Allocator>* b) {
    return *this==*b;
  }

  inline bool operator==(delegate_partitioned_graph<Allocator>& other) {
    if (m_mpi_size != other.m_mpi_size)
      return false;
    if (m_mpi_rank != other.m_mpi_rank)
      return false;
    if(m_owned_info != other.m_owned_info)
      return false;
    // if(m_owned_targets != other.m_owned_targets)
    //   return false;
    if(m_delegate_info != other.m_delegate_info)
      return false;
    if(m_delegate_degree != other.m_delegate_degree)
      return false;
    if(m_delegate_label != other.m_delegate_label)
      return false;
    // if(m_delegate_targets != other.m_delegate_targets)
    //   return false;
    if(m_map_delegate_locator != other.m_map_delegate_locator)
      return false;
    if(m_delegate_degree_threshold != other.m_delegate_degree_threshold)
      return false;
    if(m_controller_locators != other.m_controller_locators)
      return false;
    if(m_local_outgoing_count != other.m_local_outgoing_count)
      return false;
    if(m_local_incoming_count != other.m_local_incoming_count)
      return false;

    return true;
  }

  inline bool operator!=(delegate_partitioned_graph<Allocator>&
      other) {
    return !(*this == other);
  }

 private:
  //////////////////////////////////////////////////////////////////////////////
  // Private Member Functions
  //////////////////////////////////////////////////////////////////////////////
  /// Synchronizes hub set amongst all processes.
  void sync_global_hub_set(const boost::unordered_set<uint64_t>& local_hubs,
                           boost::unordered_set<uint64_t>& global_hubs,
                           bool local_change);


  void initialize_low_meta_data(boost::unordered_set<uint64_t>& global_hub_set);
  void initialize_high_meta_data(boost::unordered_set<uint64_t>& global_hubs);

  void initialize_edge_storage(const Allocator& allocator);

  template <typename Container, typename edge_data_type>
  void partition_low_degree(Container& unsorted_edges, edge_data_type& _edge_data);

  template <typename InputIterator>
  void count_high_degree_edges(InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set);


  template <typename Container, typename edge_data_type>
  void partition_high_degree(Container& unsorted_edges,
    std::map< uint64_t, std::deque<OverflowSendInfo> > &transfer_info, 
    edge_data_type& _edge_data);

  template <typename InputIterator>
  void count_edge_degrees(InputIterator unsorted_itr,
                 InputIterator unsorted_itr_end,
                 boost::unordered_set<uint64_t>& global_hub_set,
                 uint64_t delegate_degree_threshold);


  void send_high_info(std::vector< boost::container::map< uint64_t, uint64_t> >&
      maps_to_send, int maps_to_send_element_count);


  void send_vertex_info(uint64_t &high_vertex_count,
      uint64_t delegate_degree_threshold,
      std::vector<  boost::container::map<
          int, std::pair<uint64_t, uint64_t> >  >& maps_to_send,
      int maps_to_send_element_count);


  void calculate_overflow(std::map< uint64_t, std::deque<OverflowSendInfo> >
    &transfer_info);

  void generate_send_list(std::vector<uint64_t> &send_list, uint64_t num_send,
    int send_id,
    std::map< uint64_t, std::deque<OverflowSendInfo> > &transfer_info);

  void flush_graph();

  //////////////////////////////////////////////////////////////////////////////
  // Protected Data Members
  //////////////////////////////////////////////////////////////////////////////
 protected:
  uint64_t edge_chunk_size;
  int processes_per_node;
  int node_partitions;
  int m_mpi_size;
  int m_mpi_rank;
  MPI_Comm m_mpi_comm;

  ConstructionState m_graph_state {New};


  uint64_t m_max_vertex {0};
  uint64_t m_global_max_vertex {0};
  uint64_t m_global_edge_count {0};
  uint64_t m_edges_high_count {0};
  uint64_t m_edges_low_count {0};

  bip::vector<uint32_t, other_allocator<Allocator, uint32_t> > m_local_outgoing_count;
  bip::vector<uint32_t, other_allocator<Allocator, uint32_t> > m_local_incoming_count;

  bip::vector<vert_info, other_allocator<Allocator, vert_info>> m_owned_info;
  bip::vector<uint32_t, other_allocator<Allocator, uint32_t>> m_owned_info_tracker;
  //bip::vector<vertex_locator, other_allocator<Allocator, vertex_locator>> m_owned_targets;
  bip::offset_ptr<vertex_locator> m_owned_targets;
  size_t m_owned_targets_size;

  // Delegate Storage
  uint64_t m_delegate_degree_threshold;

  bip::vector< uint64_t, other_allocator<Allocator, uint64_t> > m_delegate_info;
  bip::vector< uint64_t, other_allocator<Allocator, uint64_t> > m_delegate_degree;
  bip::vector< uint64_t, other_allocator<Allocator, uint64_t> > m_delegate_label;
  // bip::vector< vertex_locator, other_allocator<Allocator, vertex_locator> >
  //     m_delegate_targets;
  bip::offset_ptr<vertex_locator> m_delegate_targets;
  size_t m_delegate_targets_size;

  //Note: BIP only contains a map, not an unordered_map object.
  /*boost::interprocess::unordered_map<
      uint64_t, vertex_locator, boost::hash<uint64_t>, std::equal_to<uint64_t>,
      other_allocator<Allocator, std::pair<uint64_t,vertex_locator> >
     > m_map_delegate_locator;
  */

  bip::map<uint64_t, vertex_locator, std::less<uint64_t>, other_allocator<Allocator, std::pair<const uint64_t,vertex_locator> > > m_map_delegate_locator;

  bip::vector<vertex_locator, other_allocator<Allocator, vertex_locator> >
    m_controller_locators;

  bool m_has_edge_data;

};  // class delegate_partitioned_graph



} // namespace havoqgt


#include <havoqgt/impl/log_step.hpp>
#include <havoqgt/impl/vert_info.hpp>
#include <havoqgt/impl/edge_data.hpp>
#include <havoqgt/impl/edge_iterator.hpp>
#include <havoqgt/impl/vertex_data.hpp>
#include <havoqgt/impl/vertex_locator.hpp>
#include <havoqgt/impl/vertex_iterator.hpp>

#include <havoqgt/impl/delegate_partitioned_graph.ipp>

#endif //HAVOQGT_MPI_DELEGATE_PARTITIONED_GRAPH_HPP_INCLUDED
