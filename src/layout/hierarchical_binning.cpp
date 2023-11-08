// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <algorithm> // for max, copy, fill_n, min
#include <cassert>   // for assert
#include <cmath>     // for log, ceil
#include <cstddef>   // for size_t
#include <limits>    // for numeric_limits
#include <string>    // for allocator, operator+, string
#include <utility>   // for pair
#include <vector>    // for vector

#include <hibf/config.hpp>                      // for config
#include <hibf/layout/data_store.hpp>           // for data_store
#include <hibf/layout/hierarchical_binning.hpp> // for hierarchical_binning
#include <hibf/layout/layout.hpp>               // for layout
#include <hibf/layout/simple_binning.hpp>       // for simple_binning
#include <hibf/misc/divide_and_ceil.hpp>        // for divide_and_ceil
#include <hibf/misc/next_multiple_of_64.hpp>    // for next_multiple_of_64
#include <hibf/platform.hpp>                    // for HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog
#include <hibf/sketch/toolbox.hpp>              // for precompute_initial_union_estimates, precompute_union_estimat...

namespace seqan::hibf::layout
{

size_t hierarchical_binning::execute()
{
    assert(data != nullptr);
    assert(data->kmer_counts != nullptr);
    assert(data->positions.size() <= data->kmer_counts->size());

    static constexpr size_t max_size_t{std::numeric_limits<size_t>::max()};

    if (!data->user_bins_arranged)
    {
        sketch::toolbox::sort_by_cardinalities(*data->kmer_counts, data->positions);

        if (!config.disable_estimate_union && !config.disable_rearrangement)
        {
            assert(data->sketches != nullptr);
            data->rearrangement_timer.start();
            sketch::toolbox::rearrange_bins(*data->sketches,
                                            *data->kmer_counts,
                                            data->positions,
                                            config.max_rearrangement_ratio,
                                            config.threads);
            data->rearrangement_timer.stop();
        }

        data->user_bins_arranged = true;
    }

    // technical bins (outer) = rows; user bins (inner) = columns
    std::vector<std::vector<size_t>> matrix(num_technical_bins, std::vector<size_t>(num_user_bins, max_size_t));

    // technical bins (outer) = rows; user bins (inner) = columns
    std::vector<std::vector<size_t>> ll_matrix(num_technical_bins, std::vector<size_t>(num_user_bins, 0u));

    // technical bins (outer) = rows; user bins (inner) = columns
    std::vector<std::vector<std::pair<size_t, size_t>>> trace(
        num_technical_bins,
        std::vector<std::pair<size_t, size_t>>(num_user_bins, {max_size_t, max_size_t}));

    initialization(matrix, ll_matrix, trace);

    recursion(matrix, ll_matrix, trace);

    // print_matrix(matrix, num_technical_bins, num_user_bins, max_size_t);
    // print_matrix(ll_matrix, num_technical_bins, num_user_bins, max_size_t);
    // print_matrix(trace, num_technical_bins, num_user_bins, std::make_pair(max_size_t, max_size_t));

    return backtracking(trace);
}

[[nodiscard]] size_t hierarchical_binning::needed_technical_bins(size_t const requested_num_ub) const
{
    return std::min<size_t>(next_multiple_of_64(requested_num_ub), config.tmax);
}

[[nodiscard]] size_t hierarchical_binning::max_merge_levels(size_t const num_ubs_in_merge) const
{
    size_t const lower_lvl_tbs = needed_technical_bins(num_ubs_in_merge);
    double const levels = std::log(num_ubs_in_merge) / std::log(lower_lvl_tbs);
    return static_cast<size_t>(std::ceil(levels));
}

void hierarchical_binning::initialization(std::vector<std::vector<size_t>> & matrix,
                                          std::vector<std::vector<size_t>> & ll_matrix,
                                          std::vector<std::vector<std::pair<size_t, size_t>>> & trace)
{
    assert(data != nullptr);

    // initialize first column
    // this is equivalent to split the first user bin into an increasing number of technical bins
    double const ub_cardinality = static_cast<double>((*data->kmer_counts)[data->positions[0]]);
    for (size_t i = 0; i < num_technical_bins; ++i)
    {
        size_t const corrected_ub_cardinality = static_cast<size_t>(ub_cardinality * data->fpr_correction[i + 1]);
        matrix[i][0] = divide_and_ceil(corrected_ub_cardinality, i + 1u);
        trace[i][0] = {0u, 0u}; // unnecessary?
    }

    // initialize first row
    // this is equivalent to merging increasing number of user bins into the first technical bin
    size_t sum = (*data->kmer_counts)[data->positions[0]];
    if (!config.disable_estimate_union)
    {
        data->union_estimation_timer.start();
        sketch::toolbox::precompute_initial_union_estimates(matrix[0],
                                                            *data->sketches,
                                                            *data->kmer_counts,
                                                            data->positions);
        data->union_estimation_timer.stop();

        for (size_t j = 1; j < num_user_bins; ++j)
        {
            sum += (*data->kmer_counts)[data->positions[j]];
            ll_matrix[0][j] = max_merge_levels(j + 1) * sum;
            trace[0][j] = {0u, j - 1}; // unnecessary?
        }
    }
    else
    {
        for (size_t j = 1; j < num_user_bins; ++j)
        {
            assert(j < data->positions.size());
            assert(data->positions[j] < data->kmer_counts->size());
            sum += (*data->kmer_counts)[data->positions[j]];
            matrix[0][j] = sum;
            ll_matrix[0][j] = max_merge_levels(j + 1) * sum;
            trace[0][j] = {0u, j - 1}; // unnecessary?
        }
    }
}

void hierarchical_binning::recursion(std::vector<std::vector<size_t>> & matrix,
                                     std::vector<std::vector<size_t>> & ll_matrix,
                                     std::vector<std::vector<std::pair<size_t, size_t>>> & trace)
{
    assert(data != nullptr);

    // we must iterate column wise
    // i iterates over the technical bins
    // j iterates over the user bins
    // matrix:
    //    i\j  UB0 UB1 UB2 UB3
    //        ---------------------------
    //    TB0 |
    //    TB1 |
    //    TB2 |
    for (size_t j = 1; j < num_user_bins; ++j)
    {
        size_t const current_weight = (*data->kmer_counts)[data->positions[j]];
        double const ub_cardinality = static_cast<double>(current_weight);

        if (!config.disable_estimate_union)
        {
            data->union_estimation_timer.start();
            sketch::toolbox::update_union_estimates_with(data->union_estimates,
                                                         *data->sketches,
                                                         *data->kmer_counts,
                                                         data->positions,
                                                         j,
                                                         config.threads);
            data->union_estimation_timer.stop();
        }

        for (size_t i = 1; i < num_technical_bins; ++i)
        {
            size_t minimum{std::numeric_limits<size_t>::max()};
            size_t full_minimum{std::numeric_limits<size_t>::max()};

            // check vertical cells
            // this is equivalent of trying to split the current user bin j into (i - i_prime) technical bins
            for (size_t i_prime = 0; i_prime < i; ++i_prime)
            {
                // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                size_t const corrected_ub_cardinality =
                    static_cast<size_t>(ub_cardinality * data->fpr_correction[(i - i_prime)]);
                size_t score =
                    std::max<size_t>(divide_and_ceil(corrected_ub_cardinality, i - i_prime), matrix[i_prime][j - 1]);
                size_t full_score = score * (i + 1) /*#TBs*/ + config.alpha * ll_matrix[i_prime][j - 1];

                // std::cout << " ++ j:" << j << " i:" << i << " i':" << i_prime << " score:" << score << std::endl;

                if (full_score < full_minimum)
                {
                    minimum = score;
                    full_minimum = full_score;
                    trace[i][j] = {i_prime, j - 1};
                    ll_matrix[i][j] = ll_matrix[i_prime][j - 1];
                }
            }

            // seqan3::debug_stream << "current vertical minimum of " << "j:" << j << " i:" << i
            //                      << " -> score:" << full_minimum << " (M_ij=" << minimum << ")"
            //                      << " trace:" << trace[i][j]
            //                      << std::endl;

            // check horizontal cells
            // this is equivalent of merging user bins {j_prime + 1, j} into technical bin i
            size_t j_prime{j - 1};
            size_t weight{current_weight};

            auto get_weight = [&]()
            {
                // if we use the union estimate we plug in that value instead of the sum (weight)
                // union_estimates[j_prime] is the union of {j_prime, ..., j}
                // the + 1 is necessary because j_prime is decremented directly after weight is updated
                return config.disable_estimate_union ? weight : data->union_estimates[j_prime + 1].estimate();
            };

            // if the user bin j-1 was not split into multiple technical bins!
            // I may merge the current user bin j into the former
            while (j_prime != 0 && ((i - trace[i][j_prime].first) < 2) && get_weight() < minimum)
            {
                weight += (*data->kmer_counts)[data->positions[j_prime]];
                --j_prime;

                // score: The current maximum technical bin size for the high-level IBF (score for the matrix M)
                // ll_kmers: estimate for the number of k-mers that have to be resolved on lower levels
                // full_score: The score to minimize -> score * #TB-high_level + low_level_memory footprint
                size_t const score = std::max<size_t>(matrix[i - 1][j_prime], get_weight());
                size_t const ll_kmers = ll_matrix[i - 1][j_prime] + max_merge_levels(j - j_prime) * weight;
                size_t const full_score = score * (i + 1) /*#TBs*/ + config.alpha * ll_kmers;

                // seqan3::debug_stream << " -- " << "j_prime:" << j_prime
                //                      << " -> full_score:" << full_score << " (M_{i-1,j'}=" << score << ")"
                //                      << std::endl;

                if (full_score < full_minimum)
                {
                    minimum = score;
                    full_minimum = full_score;
                    trace[i][j] = {i - 1, j_prime};
                    ll_matrix[i][j] = ll_kmers;
                }
            }

            matrix[i][j] = minimum;
        }
    }
}

void hierarchical_binning::backtrack_merged_bin(size_t trace_j,
                                                size_t const next_j,
                                                size_t const bin_id,
                                                maximum_bin_tracker & max_tracker,
                                                bool is_first_row)
{
    size_t kmer_count = (*data->kmer_counts)[data->positions[trace_j]];
    sketch::hyperloglog sketch =
        config.disable_estimate_union ? sketch::hyperloglog{} : (*data->sketches)[data->positions[trace_j]];
    auto libf_data = initialise_libf_data(trace_j);

    // std::cout << "merged [" << trace_j;
    trace_j -= !is_first_row;
    while (trace_j != next_j)
    {
        trace_j -= is_first_row;
        if (!config.disable_estimate_union)
            sketch.merge((*data->sketches)[data->positions[trace_j]]);
        else
            kmer_count += (*data->kmer_counts)[data->positions[trace_j]];
        libf_data.positions.push_back(data->positions[trace_j]);
        // std::cout << "," << trace_j;
        trace_j -= !is_first_row;
    }

    process_merged_bin(libf_data, bin_id);
    data->union_estimation_timer += libf_data.union_estimation_timer;
    data->rearrangement_timer += libf_data.rearrangement_timer;

    if (!config.disable_estimate_union)
        kmer_count = sketch.estimate(); // overwrite kmer_count high_level_max_id/size bin

    max_tracker.update_max(bin_id, kmer_count);
    // std::cout << "]: " << kmer_count << std::endl;
}

void hierarchical_binning::backtrack_split_bin(size_t trace_j,
                                               size_t const number_of_bins,
                                               size_t const bin_id,
                                               maximum_bin_tracker & max_tracker)
{
    assert(number_of_bins > 0u);

    // update layout
    data->hibf_layout->user_bins.emplace_back(data->previous.bin_indices,
                                              bin_id,
                                              number_of_bins,
                                              data->positions[trace_j]);

    // update max bin
    size_t const cardinality = (*data->kmer_counts)[data->positions[trace_j]];
    size_t const corrected_cardinality = static_cast<size_t>(cardinality * data->fpr_correction[number_of_bins]);
    // NOLINTNEXTLINE(clang-analyzer-core.DivideZero)
    size_t const cardinality_per_bin = divide_and_ceil(corrected_cardinality, number_of_bins);

    max_tracker.update_max(bin_id, cardinality_per_bin);
    max_tracker.update_split_max(bin_id, cardinality_per_bin);

    // std::cout << "split " << trace_j << " into " << number_of_bins << ": " << cardinality_per_bin << std::endl;
}

size_t hierarchical_binning::backtracking(std::vector<std::vector<std::pair<size_t, size_t>>> const & trace)
{
    assert(data != nullptr);

    // backtracking starts at the bottom right corner:
    size_t trace_i = num_technical_bins - 1;
    size_t trace_j = num_user_bins - 1;

    // while backtracking, keep trach of the following variables
    maximum_bin_tracker max_tracker{};
    size_t bin_id{}; // the current bin that is processed, we start naming the bins here!

    // process the trace starting at the bottom right call until you arrive at the first row or column
    while (trace_j > 0u && trace_i > 0u)
    {
        // std::cout << "\t I am now at " << trace_i << "," << trace_j << std::endl;
        size_t next_i = trace[trace_i][trace_j].first;
        size_t next_j = trace[trace_i][trace_j].second;

        size_t number_of_bins = (trace_i - next_i);

        if (number_of_bins == 1 && next_j != trace_j - 1u) // merged bin
        {
            backtrack_merged_bin(trace_j, next_j, bin_id, max_tracker);

            trace_i = next_i;
            trace_j = next_j;
        }
        else // split bin
        {
            backtrack_split_bin(trace_j, number_of_bins, bin_id, max_tracker);

            trace_i = trace[trace_i][trace_j].first;
            --trace_j;
        }

        bin_id += number_of_bins;
    }

    // process the first row or first column at last
    assert(trace_i == 0 || trace_j == 0);
    if (trace_i == 0u && trace_j > 0u) // the last UBs get merged into the remaining TB
    {
        // we are in the first row, merging the remaining UBs into the last TB (TB-0)
        backtrack_merged_bin(trace_j, 0, bin_id, max_tracker, true);
    }
    else if (trace_j == 0u) // the last UB is split into the remaining TBs
    {
        // we are in the first column, splitting the last UB (UB-0) into the remaining TBs (even if only into 1 bin).
        backtrack_split_bin(trace_j, trace_i + 1, bin_id, max_tracker);
    }

    return max_tracker.choose_max_bin(config);
}

data_store hierarchical_binning::initialise_libf_data(size_t const trace_j) const
{
    data_store libf_data{.false_positive_rate = data->false_positive_rate,
                         .hibf_layout = data->hibf_layout,
                         .kmer_counts = data->kmer_counts,
                         .sketches = data->sketches,
                         .positions = {data->positions[trace_j]},
                         .fpr_correction = data->fpr_correction};

    return libf_data;
}

void hierarchical_binning::process_merged_bin(data_store & libf_data, size_t const bin_id) const
{
    update_libf_data(libf_data, bin_id);

    // now do the binning for the low-level IBF:
    size_t const lower_max_bin = add_lower_level(libf_data);

    data->hibf_layout->max_bins.emplace_back(libf_data.previous.bin_indices, lower_max_bin);
}

void hierarchical_binning::update_libf_data(data_store & libf_data, size_t const bin_id) const
{
    bool const is_top_level = data->previous.empty();

    libf_data.previous = data->previous;
    libf_data.previous.bin_indices.push_back(bin_id);

#if HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wrestrict"
#endif // HIBF_WORKAROUND_GCC_BOGUS_MEMCPY

    libf_data.previous.num_of_bins += (is_top_level ? "" : ";") + std::string{"1"};

#if HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
#    pragma GCC diagnostic pop
#endif // HIBF_WORKAROUND_GCC_BOGUS_MEMCPY
}

size_t hierarchical_binning::add_lower_level(data_store & libf_data) const
{
    // now do the binning for the low-level IBF:
    if (libf_data.positions.size() > config.tmax)
    {
        // recursively call hierarchical binning if there are still too many UBs
        return hierarchical_binning{libf_data, config}.execute(); // return id of maximum technical bin
    }
    else
    {
        // use simple binning to distribute remaining UBs
        return simple_binning{libf_data, 0}.execute(); // return id of maximum technical bin
    }
}

} // namespace seqan::hibf::layout
