// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements seqan::hibf::compute_kmers.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cinttypes> // for uint64_t

#include <hibf/build/build_data.hpp>   // for build_data
#include <hibf/contrib/robin_hood.hpp> // for unordered_flat_set
#include <hibf/layout/layout.hpp>      // for layout

namespace seqan::hibf::build
{

void compute_kmers(robin_hood::unordered_flat_set<uint64_t> & kmers,
                   build_data const & data,
                   layout::layout::user_bin const & record);

} // namespace seqan::hibf::build
