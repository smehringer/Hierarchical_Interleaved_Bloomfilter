// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <lemon/list_graph.h> /// Must be first include.

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/detail/build/hibf/build_data.hpp>
#include <hibf/detail/build/hibf/compute_kmers.hpp>
#include <hibf/detail/build/hibf/construct_ibf.hpp>
#include <hibf/detail/build/hibf/insert_into_ibf.hpp>
#include <hibf/detail/build/hibf/loop_over_children.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <hibf/detail/build/hibf/update_user_bins.hpp>
#include <hibf/detail/build/hibf/update_parent_kmers.hpp>

namespace hibf
{

size_t hierarchical_build(robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                          lemon::ListDigraph::Node const & current_node,
                          build_data & data,
                          bool is_root)
{
    auto & current_node_data = data.node_map[current_node];

    size_t const ibf_pos{data.request_ibf_idx()};

    std::vector<int64_t> ibf_positions(current_node_data.number_of_technical_bins, ibf_pos);
    std::vector<int64_t> filename_indices(current_node_data.number_of_technical_bins, -1);
    robin_hood::unordered_flat_set<uint64_t> kmers{};

    auto initialise_max_bin_kmers = [](robin_hood::unordered_flat_set<uint64_t> & kmers,
                                       std::vector<int64_t> & ibf_positions,
                                       std::vector<int64_t> & filename_indices,
                                       lemon::ListDigraph::Node const & node,
                                       build_data & data) -> size_t
    {
        auto & node_data = data.node_map[node];

        if (node_data.favourite_child != lemon::INVALID) // max bin is a merged bin
        {
            // recursively initialize favourite child first
            ibf_positions[node_data.max_bin_index] = hierarchical_build(kmers, node_data.favourite_child, data, false);
            return 1;
        }
        else // max bin is not a merged bin
        {
            // we assume that the max record is at the beginning of the list of remaining records.
            auto const & record = node_data.remaining_records[0];
            compute_kmers(kmers, data, record);
            update_user_bins(filename_indices, record);

            return record.number_of_technical_bins;
        }
    };

    // initialize lower level IBF
    size_t const max_bin_tbs = initialise_max_bin_kmers(kmers, ibf_positions, filename_indices, current_node, data);
    auto && ibf = construct_ibf(parent_kmers, kmers, max_bin_tbs, current_node, data, is_root);
    kmers.clear(); // reduce memory peak

    // parse all other children (merged bins) of the current ibf
    auto loop_over_children = [](robin_hood::unordered_flat_set<uint64_t> & parent_kmers,
                                 hibf::interleaved_bloom_filter<> & ibf,
                                 std::vector<int64_t> & ibf_positions,
                                 lemon::ListDigraph::Node const & current_node,
                                 build_data & data,
                                 bool is_root)
    {
        auto & current_node_data = data.node_map[current_node];
        std::vector<lemon::ListDigraph::Node> children{};

        for (lemon::ListDigraph::OutArcIt arc_it(data.ibf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
            children.emplace_back(data.ibf_graph.target(arc_it));

        if (children.empty())
            return;

        size_t const number_of_mutex = (data.node_map[current_node].number_of_technical_bins + 63) / 64;
        std::vector<std::mutex> local_ibf_mutex(number_of_mutex);

        auto worker = [&](auto && index, auto &&)
        {
            auto & child = children[index];

            if (child != current_node_data.favourite_child)
            {
                robin_hood::unordered_flat_set<uint64_t> kmers{};
                size_t const ibf_pos = hierarchical_build(kmers, child, data, false);
                auto parent_bin_index = data.node_map[child].parent_bin_index;
                {
                    size_t const mutex_id{parent_bin_index / 64};
                    std::lock_guard<std::mutex> guard{local_ibf_mutex[mutex_id]};
                    ibf_positions[parent_bin_index] = ibf_pos;
                    insert_into_ibf(kmers, 1, parent_bin_index, ibf, data.fill_ibf_timer);
                    if (!is_root)
                        update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
                }
            }
        };

        size_t number_of_threads{};
        std::vector<size_t> indices(children.size());
        std::iota(indices.begin(), indices.end(), size_t{});

        if (is_root)
        {
            // Shuffle indices: More likely to not block each other. Optimal: Interleave
            std::shuffle(indices.begin(), indices.end(), std::mt19937_64{std::random_device{}()});
            number_of_threads = data.hibf_config.threads;
        }
        else
        {
            number_of_threads = 1u;
        }

        seqan3::detail::execution_handler_parallel executioner{number_of_threads};
        executioner.bulk_execute(std::move(worker), std::move(indices), []() {});
    };

    loop_over_children(parent_kmers, ibf, ibf_positions, current_node, data, is_root);

    // If max bin was a merged bin, process all remaining records, otherwise the first one has already been processed
    size_t const start{(current_node_data.favourite_child != lemon::INVALID) ? 0u : 1u};
    for (size_t i = start; i < current_node_data.remaining_records.size(); ++i)
    {
        auto const & record = current_node_data.remaining_records[i];

        if (is_root && record.number_of_technical_bins == 1) // no splitting needed
        {
            insert_into_ibf(data, record, ibf);
        }
        else
        {
            compute_kmers(kmers, data, record);
            insert_into_ibf(kmers,
                            record.number_of_technical_bins,
                            record.storage_TB_id,
                            ibf,
                            data.fill_ibf_timer);
            if (!is_root)
                update_parent_kmers(parent_kmers, kmers, data.merge_kmers_timer);
        }

        update_user_bins(filename_indices, record);
        kmers.clear();
    }

    data.hibf->ibf_vector[ibf_pos] = std::move(ibf);
    data.hibf->next_ibf_id[ibf_pos] = std::move(ibf_positions);
    data.hibf->user_bins.bin_indices_of_ibf(ibf_pos) = std::move(filename_indices);

    return ibf_pos;
}

size_t hierarchical_build(lemon::ListDigraph::Node const & root_node, build_data & data)
{
    robin_hood::unordered_flat_set<uint64_t> root_kmers{};
    return hierarchical_build(root_kmers, root_node, data, true);
}

} // namespace hibf
