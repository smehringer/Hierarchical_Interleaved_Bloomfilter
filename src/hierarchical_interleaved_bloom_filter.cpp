// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <lemon/bits/array_map.h> // for ArrayMap
#include <lemon/core.h>           // for INVALID
#include <lemon/list_graph.h>     // for ListDigraph

#include <algorithm> // for shuffle
#include <cinttypes> // for uint64_t, int64_t
#include <cstddef>   // for size_t
#include <mutex>     // for mutex, lock_guard
#include <numeric>   // for iota
#include <random>    // for random_device, mt19937_64
#include <utility>   // for move
#include <vector>    // for vector

#include <hibf/build/build_data.hpp>                      // for build_data
#include <hibf/build/compute_kmers.hpp>                   // for compute_kmers
#include <hibf/build/construct_ibf.hpp>                   // for construct_ibf
#include <hibf/build/initialise_build_tree.hpp>           // for initialise_build_tree
#include <hibf/build/insert_into_ibf.hpp>                 // for build::insert_into_ibf
#include <hibf/build/node_data.hpp>                       // for node_data
#include <hibf/build/update_parent_kmers.hpp>             // for update_parent_kmers
#include <hibf/build/update_user_bins.hpp>                // for update_user_bins
#include <hibf/config.hpp>                                // for config
#include <hibf/contrib/robin_hood.hpp>                    // for unordered_flat_set
#include <hibf/hierarchical_interleaved_bloom_filter.hpp> // for hierarchical_interleaved_bloom_filter
#include <hibf/interleaved_bloom_filter.hpp>              // for interleaved_bloom_filter
#include <hibf/layout/compute_fpr_correction.hpp>         // for compute_fpr_correction
#include <hibf/layout/compute_layout.hpp>                 // for compute_layout
#include <hibf/layout/layout.hpp>                         // for layout
#include <hibf/user_bins_type.hpp>                        // for user_bins_type

namespace seqan::hibf
{

size_t hierarchical_build(hierarchical_interleaved_bloom_filter & hibf,
                          robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                          lemon::ListDigraph::Node const & current_node,
                          build::build_data & data,
                          bool is_root)
{
    auto & current_node_data = data.node_map[current_node];

    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    robin_hood::unordered_flat_set<uint64_t> kmers{};

    auto initialise_max_bin_kmers = [&]() -> size_t
    {
        auto & node_data = data.node_map[current_node];

        if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
        {
            // recursively initialize favourite child first
            ibf_positions[node_data.max_bin_index] =
                hierarchical_build(hibf, kmers, node_data.favourite_child, data, false);
            return 1;
        }
        else // max bin is not a merged bin
        {
            // we assume that the max record is at the beginning of the list of remaining records.
            auto const & record = node_data.remaining_records[0];
            build::compute_kmers(kmers, data, record);
            build::update_user_bins(filename_indices, record);

            return record.number_of_technical_bins;
        }
    };

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers();
    auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    auto loop_over_children = [&]()
    {
        auto & current_node_data = data.node_map[current_node];
        std::vector<lemon::ListDigraph::Node> children{};

        for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            children.emplace_back(data.ibf_graph.target(arc_it));

        if (children.empty())
            return;

        size_t const number_of_mutex = (data.node_map[current_node].number_of_technical_bins + 63) / 64;
        std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

        size_t number_of_threads{};
        std::vector<size_t> indices(children.size());
        std::iota(indices.begin(), indices.end(), size_t{});

        if (is_root)
        {
            // Shuffle indices: More likely to not block each other. Optimal: Interleave
            std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
            number_of_threads = data.config.threads;
        }
        else
        {
            number_of_threads = 1u;
        }

#pragma omp parallel for schedule(dynamic, 1) num_threads(number_of_threads)
        for (auto && index : indices)
        {
            auto & child = children[index];

            if (child != current_node_data.favourite_child)
            {
                robin_hood::unordered_flat_set<uint64_t> kmers{};
                size_t const ibf_pos = hierarchical_build(hibf, kmers, child, data, false);
                auto parent_bin_index = data.node_map[child].parent_bin_index;
                {
                    size_t const mutex_id{parent_bin_index / 64};
                    std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                    ibf_positions[parent_bin_index] = ibf_pos;
                    build::insert_into_ibf(kmers, 1, parent_bin_index, ibf, data.fill_ibf_timer);
                    if (!is_root)
                        build::update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
                }
            }
        }
    };

    loop_over_children();

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];

        if (is_root && record.number_of_technical_bins == 1) // no splitting needed
        {
            build::insert_into_ibf(data, record, ibf);
        }
        else
        {
            compute_kmers(kmers, data, record);
            build::insert_into_ibf(kmers, record.number_of_technical_bins, record.storage_TB_id, ibf, data.fill_ibf_timer);
            if (!is_root)
                build::update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
        }

        build::update_user_bins(filename_indices, record);
        kmers.clear();
    }

    hibf.ibf_vector[ibf_pos] = std::move(ibf);
    hibf.next_ibf_id[ibf_pos] = std::move(ibf_positions);
    hibf.user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

    return ibf_pos;
}

size_t hierarchical_build(hierarchical_interleaved_bloom_filter & hibf,
                          lemon::ListDigraph::Node const & root_node,
                          build::build_data & data)
{
    robin_hood::unordered_flat_set<uint64_t> root_kmers{};
    return hierarchical_build(hibf, root_kmers, root_node, data, true);
}

void build_index(hierarchical_interleaved_bloom_filter & hibf,
                 config const & config,
                 seqan::hibf::layout::layout & hibf_layout)
{
    size_t const number_of_ibfs = hibf_layout.max_bins.size() + 1;

    hibf.ibf_vector.resize(number_of_ibfs);
    hibf.user_bins.set_ibf_count(number_of_ibfs);
    hibf.user_bins.set_user_bin_count(hibf_layout.user_bins.size());
    hibf.next_ibf_id.resize(number_of_ibfs);

    build::build_data data{.config = config};

    initialise_build_tree(hibf_layout, data.ibf_graph, data.node_map);
    lemon::ListDigraph::Node root_node = data.ibf_graph.nodeFromId(0); // root node = top-level IBF node

    size_t const t_max{data.node_map[root_node].number_of_technical_bins};
    data.fpr_correction = layout::compute_fpr_correction(
        {.fpr = config.maximum_false_positive_rate, .hash_count = config.number_of_hash_functions, .t_max = t_max});

    hierarchical_build(hibf, root_node, data);
}

hierarchical_interleaved_bloom_filter::hierarchical_interleaved_bloom_filter(config const & configuration)
{
    auto layout = layout::compute_layout(configuration);
    build_index(*this, configuration, layout);
}

hierarchical_interleaved_bloom_filter::hierarchical_interleaved_bloom_filter(
    std::function<void(size_t const, insert_iterator &&)> input_fn,
    std::istream & layout_stream)
{
    // read config and layout from file
    config configuration;
    layout::layout hibf_layout;
    configuration.read_from(layout_stream);
    hibf_layout.read_from(layout_stream);

    configuration.input_fn = input_fn; // set input as it cannot be serialized.

    build_index(*this, configuration, hibf_layout);
}

} // namespace seqan::hibf
