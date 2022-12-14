
// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#ifndef HAVOQGT_MPI_IMPL_VERTEX_LOCATOR_HPP_
#define HAVOQGT_MPI_IMPL_VERTEX_LOCATOR_HPP_

#include <havoqgt/delegate_partitioned_graph.hpp>

namespace havoqgt {

template <typename SegementManager>
class delegate_partitioned_graph<SegementManager>::vertex_locator {
 public:
  vertex_locator() {
    m_is_delegate  = 0;
    m_is_bcast     = 0;
    m_is_intercept = 0;
#pragma GCC diagnostic ignored \
    "-Woverflow"  /// NOTE:  is there a better way to clean these overflows?
    m_owner_dest = std::numeric_limits<uint32_t>::max();
    m_local_id   = std::numeric_limits<uint64_t>::max();
#pragma GCC diagnostic pop
  }

  bool is_valid() const {
    delegate_partitioned_graph<SegementManager>::vertex_locator conv;
#pragma GCC diagnostic ignored \
    "-Woverflow"  /// NOTE:  is there a better way to clean these overflows?
    conv.m_local_id   = std::numeric_limits<uint64_t>::max();
    conv.m_owner_dest = std::numeric_limits<uint64_t>::max();
#pragma GCC diagnostic pop

    return (m_local_id != conv.m_local_id || m_owner_dest != conv.m_owner_dest);
  }
  bool is_delegate_master() const {
    return (is_delegate() &&
            ((m_local_id % comm_world().size()) == comm_world().rank()));
  }

  bool     is_delegate() const { return m_is_delegate == 1; }
  uint32_t owner() const { return m_owner_dest; }
  void set_dest(uint32_t dest) {
    m_owner_dest = dest;
    assert(m_owner_dest == dest);
  }
  uint64_t local_id() const { return m_local_id; }
  bool is_equal(const vertex_locator x) const;
  uint32_t get_bcast() const { return m_is_bcast; }
  void set_bcast(uint32_t bcast) { m_is_bcast = bcast; }
  bool                    is_intercept() const { return m_is_intercept == 1; }
  void set_intercept(bool intercept) { m_is_intercept = intercept; }

  uint64_t hash() const {
    // return (uint64_t(m_owner_dest) << 39) | uint64_t(m_local_id) |
    //        (uint64_t(m_is_delegate) << 60);
    uint64_t key;
    if (m_is_delegate) {
      key = (uint64_t(1) << 62) | uint64_t(m_local_id);
    } else {
      key = (uint64_t(m_owner_dest) << 40) | uint64_t(m_local_id);
    }
    return hash64shift(key);
  }

  friend bool operator==(const vertex_locator& x, const vertex_locator& y) {
    return x.is_equal(y);
  }
  friend bool operator<(const vertex_locator& x, const vertex_locator& y) {
    //  THIS WIS BUGGY WITH delegate to delegate comparisons
    // if (x.m_is_delegate &&
    //     x.m_owner_dest != x.m_local_id % comm_world().size()) {
    //   std::cerr << "X owner_dest for delegate issue" << std::endl;
    //   exit(-1);
    // }
    // if (y.m_is_delegate &&
    //     y.m_owner_dest != y.m_local_id % comm_world().size()) {
    //   std::cerr << "Y owner_dest for delegate issue" << std::endl;
    //   exit(-1);
    // }
    // if (x.m_is_delegate == y.m_is_delegate) {
    //   if (x.m_owner_dest == y.m_owner_dest) {
    //     return x.m_local_id < y.m_local_id;
    //   } else {
    //     return x.m_owner_dest < y.m_owner_dest;
    //   }

    // } else {
    //   return x.m_is_delegate < y.m_is_delegate;
    // }

    if (x.m_is_delegate && y.m_is_delegate) {
      // both delegates
      return x.m_local_id < y.m_local_id;
    } else if (!x.m_is_delegate && !y.m_is_delegate) {
      // both not delegates
      if (x.m_owner_dest < y.m_owner_dest) {
        return true;
      } else if (x.m_owner_dest > y.m_owner_dest) {
        return false;
      } else {
        return x.m_local_id < y.m_local_id;
      }
    } else {
      return x.m_is_delegate < y.m_is_delegate;
    }
  }

  friend bool operator!=(const vertex_locator& x, const vertex_locator& y) {
    return !(x.is_equal(y));
  }

 private:
  uint64_t hash64shift(uint64_t key) const {
    key = (~key) + (key << 21);  // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);  // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);  // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
  }

  friend class delegate_partitioned_graph;
  unsigned int m_is_delegate : 1;

  unsigned int m_is_bcast : 1;
  unsigned int m_is_intercept : 1;
  unsigned int m_owner_dest : 21;
  uint64_t     m_local_id : 40;

  vertex_locator(bool is_delegate, uint64_t local_id, uint32_t owner_dest);
} __attribute__((packed));

///////////////////////////////////////////////////////////////////////////////
//                           Vertex Locator                                  //
///////////////////////////////////////////////////////////////////////////////
/**
 * @class  delegate_partitioned_graph::vertex_locator
 * @details Here are some very important details.
 */
/**
 *
 */
template <typename SegmentManager>
inline delegate_partitioned_graph<
    SegmentManager>::vertex_locator::vertex_locator(bool     is_delegate,
                                                    uint64_t local_id,
                                                    uint32_t owner_dest) {
  m_is_bcast     = 0;
  m_is_intercept = 0;

  if (is_delegate) {
    m_is_delegate = true;
    m_owner_dest  = owner_dest;
    m_local_id    = local_id;
    if (!(m_is_delegate == true && m_local_id == local_id &&
          m_owner_dest == owner_dest)) {
      std::cerr << "ERROR:  vertex_locator()" << std::endl;
      exit(-1);
    }
  } else {
    m_is_delegate = false;
    m_owner_dest  = owner_dest;
    m_local_id    = local_id;
    if (!(m_is_delegate == false && m_owner_dest == owner_dest &&
          m_local_id == local_id)) {
      std::cerr << "ERROR:  vertex_locator()" << std::endl;
      exit(-1);
    }
  }
}

template <typename SegmentManager>
inline bool
delegate_partitioned_graph<SegmentManager>::vertex_locator::is_equal(
    const typename delegate_partitioned_graph<SegmentManager>::vertex_locator x)
    const {
  return m_is_delegate == x.m_is_delegate
         //&& m_is_bcast == x.m_is_bcast &&
         // m_is_intercept == x.m_is_intercept
         && m_owner_dest == x.m_owner_dest && m_local_id == x.m_local_id;
}

}  // namespace havoqgt

#endif  // HAVOQGT_MPI_IMPL_VERTEX_LOCATOR_HPP_
