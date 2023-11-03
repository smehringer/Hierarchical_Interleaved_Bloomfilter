// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cassert> // for assert
#include <cstddef> // for size_t
#include <utility> // for addressof, pair
#include <vector>  // for vector

#include <hibf/build/bin_size_in_bits.hpp> // for bin_size_in_bits
#include <hibf/config.hpp>                 // for config
#include <hibf/layout/data_store.hpp>      // for data_store
#include <hibf/sketch/hyperloglog.hpp>     // for hyperloglog
#include <hibf/platform.hpp>               // for HIBF_WORKAROUND_GCC_BOGUS_MEMCPY

namespace seqan::hibf::layout
{

/*!\brief Hierarchical binning algorithm.
 * \ingroup hibf_layout
 */
class hierarchical_binning
{
private:
    //!\brief The user config passed down from the command line.
    seqan::hibf::config config{};
    //!\brief Stores all data that is needed to compute the layout, e.g. the counts, sketches and the layout::layout.
    data_store * data{nullptr};

    //!\brief The number of user bins, initialised with the length of user_bin_kmer_counts.
    size_t num_user_bins{};
    //!\brief The number of technical bins requested by the user.
    size_t num_technical_bins{};

    /*!\brief given the number of unique kmers in the input and the number of TBs to distribute, in an optimal layout,
     *        a technical bin will never have a lower weight than `minimum_technical_bin_weight`.
     */
    size_t minimum_technical_bin_weight{0};

    //!\brief Simplifies passing the parameters needed for tracking the maximum technical bin.
    struct maximum_bin_tracker
    {
        size_t max_id{};         //!< The ID of the technical bin with maximal size.
        size_t max_size{};       //!< The maximum technical bin size seen so far.
        size_t max_split_id{};   //!< The ID of the split bin with maximal size (if any).
        size_t max_split_size{}; //!< The maximum split bin size seen so far.

        void update_max(size_t const new_id, size_t const new_size)
        {
            if (new_size > max_size)
            {
                max_id = new_id;
                max_size = new_size;
            }
        }

        //!\brief Split cardinality `new_size` must already account for fpr-correction.
        void update_split_max(size_t const new_id, size_t const new_size)
        {
            if (new_size > max_split_size)
            {
                max_split_id = new_id;
                max_split_size = new_size;
            }
        }

        /*!\brief Decides which bin is reported as the maximum bin.
         *\param config The HIBF configuration.
         *\return The chosen max bin id.
         *
         * As a HIBF feature, the merged bin FPR can differ from the overall maximum FPR. Merged bins in an HIBF layout
         * will always be followed by one or more lower-level IBFs that will have split bins or single bins (split = 1)
         * to recover the original user bins.
         *
         * We need to make sure, though, that downsizing merged bins does not affect split bins.
         * Therefore, we check if choosing a merged bin as the max bin violates the minimum_bits needed for split bins.
         * If so, we can report the largest split bin as the max bin as it will choose the correct size and downsize
         * larger merged bins only a little.
         */
        size_t choose_max_bin(seqan::hibf::config const & config)
        {
            if (max_id == max_split_id) // Overall max bin is a split bin.
                return max_id;

            // Split cardinality `max_split_size` already accounts for fpr correction.
            // The minimum size of the TBs of this IBF to ensure the maximum_false_positive_rate for split bins.
            size_t const minimum_bits{build::bin_size_in_bits({.fpr = config.maximum_fpr,
                                                               .hash_count = config.number_of_hash_functions,
                                                               .elements = max_split_size})};

            // The potential size of the TBs of this IBF given the allowed merged bin FPR.
            size_t const merged_bits{build::bin_size_in_bits({.fpr = config.relaxed_fpr, //
                                                              .hash_count = config.number_of_hash_functions,
                                                              .elements = max_size})};

            // If split and merged bits are the same, we prefer merged bins. Better for build parallelisation.
            if ((minimum_bits > merged_bits))
                return max_split_id;

            return max_id;
        }
    };

public:
    hierarchical_binning() = default;                                        //!< Defaulted.
    hierarchical_binning(hierarchical_binning const &) = delete;             //!< Deleted. Would modify same data.
    hierarchical_binning & operator=(hierarchical_binning const &) = delete; //!< Deleted. Would modify same data.
    hierarchical_binning(hierarchical_binning &&) = default;                 //!< Defaulted.
    hierarchical_binning & operator=(hierarchical_binning &&) = default;     //!< Defaulted.
    ~hierarchical_binning() = default;                                       //!< Defaulted.

    /*!\brief The constructor from user bin names, their kmer counts and a config.
     * \param[in, out] data_ Stores all data that is needed to compute the layout.
     * \param[in] config_ A config object that holds information from the user that influence the computation.
     *
     * Each entry in the names_ and input vector respectively is considered a user bin (both vectors must have the
     * same length).
     */
    hierarchical_binning(data_store & data_, seqan::hibf::config const & config_) :
        config{config_},
        data{std::addressof(data_)},
        num_user_bins{data->positions.size()},
        num_technical_bins{data->previous.empty() ? config.tmax : needed_technical_bins(num_user_bins)}
    {
        assert(data != nullptr);
        sketch::hyperloglog sketch(config.sketch_bits);
        assert((data->sketches));
        for (auto const pos : data->positions)
            sketch.merge((*data->sketches)[pos]);
        minimum_technical_bin_weight = std::max(sketch.estimate() / num_technical_bins, 1.0);
    }

    //!\brief Executes the hierarchical binning algorithm and layouts user bins into technical bins.
    size_t execute();

private:
    /*!\brief Returns the number of technical bins given a number of user bins.
     * \param[in] requested_num_ub The number of user bins.
     */
    [[nodiscard]] size_t needed_technical_bins(size_t const requested_num_ub) const;

    /*!\brief Returns the maximum number of needed levels when merging `num_ubs_in_merge` many user bins.
     * \param[in] num_ubs_in_merge The number of user bins in the merge.
     */
    [[nodiscard]] size_t max_merge_levels(size_t const num_ubs_in_merge) const;

    /*!\brief Initialize the matrices M (hig_level_ibf), L (low_level_ibfs) and T (trace)
     *
     * \image html hierarchical_dp_init.png
     */
    void initialization(std::vector<std::vector<size_t>> & matrix,
                        std::vector<std::vector<size_t>> & ll_matrix,
                        std::vector<std::vector<std::pair<size_t, size_t>>> & trace);

    /*!\brief Performs the recursion.
     *
     * \image html hierarchical_dp_recursion.png
     *
     * Explanations to the formula:
     *
     * Remember that M (matrix) stores the maximum technical bin size of the high level IBF (HIBF) and
     * L (ll_matrix) stores the sum of all estimate low level ibfs (LIBF) memory footprints
     * (assuming perfect splitting of kmer_content which can obviously not be achieved but `alpha` may be adjusted
     * in order to counteract this underestimation).
     *
     * Now in order to minimize the memory footprint we...
     *
     * 1. ... (\f$v_ij\f$) firstly check if we should split the bin. <br>
     * Therefore we compute for every possible \f$ i' \f$ above \f$ i \f$, the technical bin size if we split
     * \f$ c_j \f$ into \f$ i - i' \f$ bins (\f$ \frac{c_j}{i - i'} \f$). We only take the maximum of the current
     * maximum where I come from (\f$ M_{i',j-1} \f$) and the new technical bin size computed just now.
     * This maximum is the new current *maximum technical bin size* that needs to be multiplied by the current number
     * technical bins \f$ (i + 1) \f$ in order to estimate the memory footprint of the HIBF.
     * Now that we have the memory footprint of the HIBF we also consider the LIBFs memory footprint but nothing
     * changed here, since we split not merge, so we just take \f$ L_{i',j-1} \f$ scaled by alpha.
     *
     * 2. ... (\f$h_ij\f$) secondly check if we should merge the bin with the ones before.<br>
     * Therefore we start by compute for every possible \f$ j' \f$ to the left of \f$ j \f$, the merged bin weight
     * (\f$ \sum_{g = j'}^{j} c_g \f$) of merging together user bins \f$ [j', ...,  j] \f$. This is only possible
     * iff every bin \f$ [j', ...,  j] \f$ was **not splitted**. If we merge those bins starting with user bin
     * \f$ j' \f$ we start the trace at \f$ M_{i-1,j'-1} \f$. Therefore we need to compute the new maximal technical
     * bin size of the HIBF, which is the maximum of the merged bin weight and where we would come from
     * (\f$ \max(M_{i-1,j'-1}, \sum_{g = j'}^{j} c_g) \f$) and multiple it by the number of technical bins so far
     * \f$ (i + 1) \f$ to get the HIBF memory footprint. The LIBFs memory footprint also changes since we introduce a
     * new merged bin. Namely, we add the weight of the new merged bin (\f$ \sum_{g = j'}^{j} c_g \f$) again to the
     * LIBFs memory footprint from where we would come from \f$ L_{i-1,j'-1} \f$ and scale this by alpha. Just adding
     * the combined merged bin weight neglects the fact, that the merged bin weight has to be distributed within the
     * new low level IBF, potentially causing the effective text ratio and thereby the IBF memory footprint to increase.
     * This we cannot know before hand how the data are, we need to accept this as a flaw in the "optimumal result" of
     * this algorithm. It would be too computational intensive to compute the splitting for every possibility.
     *
     */
    void recursion(std::vector<std::vector<size_t>> & matrix,
                   std::vector<std::vector<size_t>> & ll_matrix,
                   std::vector<std::vector<std::pair<size_t, size_t>>> & trace);

    void backtrack_merged_bin(size_t trace_j,
                              size_t const next_j,
                              size_t const bin_id,
                              maximum_bin_tracker & max_tracker,
                              bool is_first_row = false);

    void backtrack_split_bin(size_t trace_j,
                             size_t const number_of_bins,
                             size_t const bin_id,
                             maximum_bin_tracker & max_tracker);

    //!\brief Backtracks the trace matrix and writes the resulting binning into the output file.
    size_t backtracking(std::vector<std::vector<std::pair<size_t, size_t>>> const & trace);

    data_store initialise_libf_data(size_t const trace_j) const;

    void process_merged_bin(data_store & libf_data, size_t const bin_id) const;

    void update_libf_data(data_store & libf_data, size_t const bin_id) const;

    size_t add_lower_level(data_store & libf_data) const;
};

} // namespace seqan::hibf::layout
