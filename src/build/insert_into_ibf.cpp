// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <cassert>    // for assert
#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iterator>   // for inserter, counted_iterator
#include <ranges>     // for all_t, operator==

#include <hibf/build/build_data.hpp>                // for build_data
#include <hibf/build/insert_into_ibf.hpp>           // for insert_into_ibf
#include <hibf/config.hpp>                          // for config
#include <hibf/contrib/robin_hood.hpp>              // for unordered_flat_set
#include <hibf/contrib/std/chunk_view.hpp>          // for chunk_view, operator==, chunk, chunk_fn
#include <hibf/contrib/std/detail/adaptor_base.hpp> // for operator|
#include <hibf/detail/timer.hpp>                    // for concurrent, timer
#include <hibf/interleaved_bloom_filter.hpp>        // for interleaved_bloom_filter, bin_index
#include <hibf/layout/layout.hpp>                   // for layout

namespace seqan::hibf::build
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan::hibf::interleaved_bloom_filter & ibf,
                     timer<concurrent::yes> & fill_ibf_timer)
{
    size_t const chunk_size = kmers.size() / number_of_bins + 1;
    size_t chunk_number{};

    timer<concurrent::no> local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    for (auto chunk : kmers | seqan::stl::views::chunk(chunk_size))
    {
        assert(chunk_number < number_of_bins);
        seqan::hibf::bin_index const bin_idx{bin_index + chunk_number};
        ++chunk_number;
        for (size_t const value : chunk)
            ibf.emplace(value, bin_idx);
    }
    local_fill_ibf_timer.stop();
    fill_ibf_timer += local_fill_ibf_timer;
}

void insert_into_ibf(build_data const & data,
                     layout::layout::user_bin const & record,
                     seqan::hibf::interleaved_bloom_filter & ibf)
{
    auto const bin_index = seqan::hibf::bin_index{static_cast<size_t>(record.storage_TB_id)};
    robin_hood::unordered_flat_set<uint64_t> values;

    timer<concurrent::no> local_user_bin_io_timer{};
    local_user_bin_io_timer.start();
    data.config.input_fn(record.idx, std::inserter(values, values.begin()));
    local_user_bin_io_timer.stop();
    data.user_bin_io_timer += local_user_bin_io_timer;

    timer<concurrent::no> local_fill_ibf_timer{};
    local_fill_ibf_timer.start();
    for (auto && value : values)
        ibf.emplace(value, bin_index);
    local_fill_ibf_timer.stop();
    data.fill_ibf_timer += local_fill_ibf_timer;
}

} // namespace seqan::hibf::build
