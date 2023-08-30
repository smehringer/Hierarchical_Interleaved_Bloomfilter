// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <lemon/list_graph.h> // for ListDigraph

#include <cinttypes> // for uint64_t
#include <cstddef>   // for size_t

#include <hibf/build/build_data.hpp>         // for build_data
#include <hibf/contrib/robin_hood.hpp>       // for unordered_flat_set
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter

namespace seqan::hibf::build
{

seqan::hibf::interleaved_bloom_filter construct_ibf(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                                    robin_hood::unordered_flat_set<uint64_t> & kmers,
                                                    size_t const number_of_bins,
                                                    lemon::ListDigraph::Node const & node,
                                                    build_data & data,
                                                    bool is_root);

} // namespace seqan::hibf::build
