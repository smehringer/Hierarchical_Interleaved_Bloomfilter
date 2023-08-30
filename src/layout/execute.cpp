#include <cinttypes> // for uint16_t
#include <cmath>     // for ceil, sqrt
#include <cstddef>   // for size_t
#include <iostream>  // for char_traits, operator<<, basic_ostream, cerr
#include <stdexcept> // for invalid_argument
#include <vector>    // for vector

#include <hibf/config.hpp>                        // for config
#include <hibf/layout/compute_fpr_correction.hpp> // for compute_fpr_correction
#include <hibf/layout/data_store.hpp>             // for data_store
#include <hibf/layout/execute.hpp>                // for execute
#include <hibf/layout/hierarchical_binning.hpp>   // for hierarchical_binning
#include <hibf/next_multiple_of_64.hpp>           // for next_multiple_of_64

namespace seqan::hibf::layout
{

size_t execute(seqan::hibf::config const & config, seqan::hibf::layout::data_store & data)
{
    seqan::hibf::config config_copy{config};

    if (config_copy.disable_estimate_union)
        config_copy.disable_rearrangement = true;

    if (config_copy.tmax == 0) // no tmax was set by the user on the command line
    {
        // Set default as sqrt(#samples). Experiments showed that this is a reasonable default.
        if (size_t number_samples = data.kmer_counts->size();
            number_samples >= 1ULL << 32) // sqrt is bigger than uint16_t
            throw std::invalid_argument{"Too many samples. Please set a tmax (see help via `-hh`)."}; // GCOVR_EXCL_LINE
        else
            config_copy.tmax =
                seqan::hibf::next_multiple_of_64(static_cast<uint16_t>(std::ceil(std::sqrt(number_samples))));
    }
    else if (config_copy.tmax % 64 != 0)
    {
        config_copy.tmax = seqan::hibf::next_multiple_of_64(config_copy.tmax);
        std::cerr << "[HIBF LAYOUT WARNING]: Your requested number of technical bins was not a multiple of 64. "
                  << "Due to the architecture of the HIBF, it will use up space equal to the next multiple of 64 "
                  << "anyway, so we increased your number of technical bins to " << config_copy.tmax << ".\n";
    }

    data.fpr_correction =
        compute_fpr_correction({.fpr = config_copy.maximum_false_positive_rate, // prevent clang-format
                                .hash_count = config_copy.number_of_hash_functions,
                                .t_max = config_copy.tmax});

    return seqan::hibf::layout::hierarchical_binning{data, config_copy}.execute();
}

} // namespace seqan::hibf::layout
