// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes> // for uint64_t
#include <cstddef>   // for size_t

#include <hibf/build/build_data.hpp>         // for build_data
#include <hibf/contrib/robin_hood.hpp>       // for unordered_flat_set
#include <hibf/detail/timer.hpp>             // for concurrent, timer
#include <hibf/interleaved_bloom_filter.hpp> // for interleaved_bloom_filter
#include <hibf/layout/layout.hpp>            // for layout

namespace seqan::hibf::build
{

// automatically does naive splitting if number_of_bins > 1
void insert_into_ibf(robin_hood::unordered_flat_set<uint64_t> const & kmers,
                     size_t const number_of_bins,
                     size_t const bin_index,
                     seqan::hibf::interleaved_bloom_filter & ibf,
                     timer<concurrent::yes> & fill_ibf_timer);

void insert_into_ibf(build_data const & data,
                     layout::layout::user_bin const & record,
                     seqan::hibf::interleaved_bloom_filter & ibf);

} // namespace seqan::hibf::build
