// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

#include <unordered_map>
#include "config.hpp"

namespace prunejuice {
#ifdef ENABLE_MY_LCC
template<typename Vertex, typename VertexData, typename BitSet, typename NLFMap>
#else
template<typename Vertex, typename VertexData, typename BitSet>
#endif
class vertex_state {
  public :
    vertex_state() :
    vertex_pattern_index(0),
    is_active(false) {}
  
    BitSet template_vertices;
    BitSet template_neighbors;
#ifdef ENABLE_MY_LCC
    NLFMap neighbor_label_frequency_map;
#endif
    //std::unordered_map<VertexData, IntegralType>
    //  template_neighbor_metadata_count_map;
 
    size_t vertex_pattern_index; // TODO: dummy, to be removed
    bool is_active; 		  
};
  
} // end namespace prunejuice 
