// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/bits/array_map.h> // for ArrayMap
#include <lemon/list_graph.h>     // for ListDigraph

#include <cinttypes> // for uint64_t
#include <cmath>     // for ceil
#include <cstddef>   // for size_t
#include <vector>    // for vector

#include <hibf/build/bin_size_in_bits.hpp>    // for bin_size_in_bits
#include <hibf/build/build_data.hpp>          // for build_data
#include <hibf/build/construct_ibf.hpp>       // for construct_ibf
#include <hibf/build/insert_into_ibf.hpp>     // for insert_into_ibf
#include <hibf/build/node_data.hpp>           // for node_data
#include <hibf/build/update_parent_kmers.hpp> // for update_parent_kmers
#include <hibf/config.hpp>                    // for config
#include <hibf/contrib/robin_hood.hpp>        // for unordered_flat_set
#include <hibf/detail/timer.hpp>              // for concurrent, timer
#include <hibf/interleaved_bloom_filter.hpp>  // for interleaved_bloom_filter, bin_count, bin_size, hash_fun...

namespace seqan::hibf::build
{

seqan::hibf::interleaved_bloom_filter construct_ibf(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                                    robin_hood::unordered_flat_set<uint64_t> & kmers,
                                                    size_t const number_of_bins,
                                                    lemon::ListDigraph::Node const & node,
                                                    build_data & data,
                                                    bool is_root)
{
    auto & node_data = data.node_map[node];

    size_t const kmers_per_bin{static_cast<size_t>(std::ceil(static_cast<double>(kmers.size()) / number_of_bins))};
    double const bin_bits{static_cast<double>(bin_size_in_bits({.fpr = data.config.maximum_false_positive_rate,
                                                                .hash_count = data.config.number_of_hash_functions,
                                                                .elements = kmers_per_bin}))};
    seqan::hibf::bin_size const bin_size{
        static_cast<size_t>(std::ceil(bin_bits * data.fpr_correction[number_of_bins]))};
    seqan::hibf::bin_count const bin_count{node_data.number_of_technical_bins};

    timer<concurrent::no> local_index_allocation_timer{};
    local_index_allocation_timer.start();
    seqan::hibf::interleaved_bloom_filter ibf{bin_count,
                                              bin_size,
                                              seqan::hibf::hash_function_count{data.config.number_of_hash_functions}};
    local_index_allocation_timer.stop();
    data.index_allocation_timer += local_index_allocation_timer;

    insert_into_ibf(kmers, number_of_bins, node_data.max_bin_index, ibf, data.fill_ibf_timer);
    if (!is_root)
        update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);

    return ibf;
}

} // namespace seqan::hibf::build
