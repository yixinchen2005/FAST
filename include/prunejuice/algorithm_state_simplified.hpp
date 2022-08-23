// Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
// HavoqGT Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: MIT

#pragma once

namespace prunejuice {

template<typename BitSet, typename NLFMap>
class vertex_state {
  public :
    vertex_state() {}
  
    BitSet template_candidates;
    NLFMap neighbor_label_frequency_map;		  
};
  
} // end namespace prunejuice 