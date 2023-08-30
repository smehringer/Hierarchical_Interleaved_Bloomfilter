#include <cinttypes>  // for uint64_t
#include <cstddef>    // for size_t
#include <functional> // for function
#include <iosfwd>     // for stringstream
#include <iterator>   // for inserter
#include <utility>    // for addressof
#include <vector>     // for vector

#include <hibf/config.hpp>                      // for config
#include <hibf/contrib/robin_hood.hpp>          // for unordered_flat_set
#include <hibf/layout/compute_layout.hpp>       // for compute_layout
#include <hibf/layout/data_store.hpp>           // for data_store
#include <hibf/layout/execute.hpp>              // for execute
#include <hibf/layout/layout.hpp>               // for layout
#include <hibf/sketch/estimate_kmer_counts.hpp> // for estimate_kmer_counts
#include <hibf/sketch/hyperloglog.hpp>          // for hyperloglog

namespace seqan::hibf::layout
{

layout
compute_layout(config const & config, std::vector<size_t> & kmer_counts, std::vector<sketch::hyperloglog> & sketches)
{
    layout resulting_layout{};

    // The output streams facilitate writing the layout file in hierarchical structure.
    // seqan::hibf::execute currently writes the filled buffers to the output file.
    std::stringstream output_buffer;
    std::stringstream header_buffer;

    // compute sketches
    sketches.resize(config.number_of_user_bins);
    kmer_counts.resize(config.number_of_user_bins);

    robin_hood::unordered_flat_set<uint64_t> kmers;
#pragma omp parallel for schedule(static) num_threads(config.threads) private(kmers)
    for (size_t i = 0; i < config.number_of_user_bins; ++i)
    {
        seqan::hibf::sketch::hyperloglog sketch(config.sketch_bits);

        kmers.clear();
        config.input_fn(i, std::inserter(kmers, kmers.begin()));

        for (auto k_hash : kmers)
            sketch.add(reinterpret_cast<char *>(&k_hash), sizeof(k_hash));

        // #pragma omp critical
        sketches[i] = sketch;
    }

    sketch::estimate_kmer_counts(sketches, kmer_counts);

    data_store store{.false_positive_rate = config.maximum_false_positive_rate,
                     .hibf_layout = &resulting_layout,
                     .kmer_counts = std::addressof(kmer_counts),
                     .sketches = std::addressof(sketches)};

    size_t const max_hibf_id = seqan::hibf::layout::execute(config, store);
    store.hibf_layout->top_level_max_bin_id = max_hibf_id;

    return *store.hibf_layout; // return layout as string for now, containing the file
}

layout compute_layout(config const & config)
{
    std::vector<size_t> kmer_counts{};
    std::vector<sketch::hyperloglog> sketches{};

    return compute_layout(config, kmer_counts, sketches);
}

} // namespace seqan::hibf::layout
