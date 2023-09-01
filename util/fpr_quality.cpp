// ------------------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/Hierarchical_Interleaved_Bloomfilter/blob/main/LICENSE.md
// ------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <random>

#include <hibf/build/bin_size_in_bits.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/interleaved_bloom_filter.hpp>
#include <hibf/layout/compute_fpr_correction.hpp>

#include <sharg/parser.hpp>

struct config
{
    // Set by user
    size_t kmer_size{12};
    size_t elements{130};
    size_t splits{5};
    size_t hash{2};
    double fpr{0.05};

    // Internal
    size_t number_of_kmers{};
    size_t split_elements_per_bin{};
};

class positive_integer_validator
{
public:
    using option_value_type = size_t;

    positive_integer_validator() = default;
    positive_integer_validator(positive_integer_validator const &) = default;
    positive_integer_validator & operator=(positive_integer_validator const &) = default;
    positive_integer_validator(positive_integer_validator &&) = default;
    positive_integer_validator & operator=(positive_integer_validator &&) = default;
    ~positive_integer_validator() = default;

    explicit positive_integer_validator(bool const is_zero_positive_) : is_zero_positive{is_zero_positive_}
    {}

    void operator()(option_value_type const & val) const
    {
        if (!is_zero_positive && !val)
            throw sharg::validation_error{"The value must be a positive integer."};
    }

    std::string get_help_page_message() const
    {
        if (is_zero_positive)
            return "Value must be a positive integer or 0.";
        else
            return "Value must be a positive integer.";
    }

private:
    bool is_zero_positive{false};
};

void init_parser(sharg::parser & parser, config & cfg)
{
    parser.add_option(cfg.kmer_size,
                      sharg::config{.short_id = '\0',
                                    .long_id = "kmer",
                                    .description = "The k-mer size.",
                                    .validator = sharg::arithmetic_range_validator{1, 31}});
    parser.add_option(cfg.elements,
                      sharg::config{.short_id = '\0',
                                    .long_id = "elements",
                                    .description = "Number of elements to insert.",
                                    .validator = positive_integer_validator{}});
    parser.add_option(cfg.splits,
                      sharg::config{.short_id = '\0',
                                    .long_id = "splits",
                                    .description = "Number of bins to split into.",
                                    .validator = positive_integer_validator{}});
    parser.add_option(cfg.hash,
                      sharg::config{.short_id = '\0',
                                    .long_id = "hash",
                                    .description = "The number of hash functions to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 5}});
    parser.add_option(cfg.fpr,
                      sharg::config{.short_id = '\0',
                                    .long_id = "fpr",
                                    .description = "The desired false positive rate.",
                                    .validator = sharg::arithmetic_range_validator{0.0, 1.0}});
}

size_t split_bin_size_in_bits(config const & cfg)
{
    return seqan::hibf::build::bin_size_in_bits(
        {.fpr = cfg.fpr, .hash_count = cfg.hash, .elements = cfg.split_elements_per_bin});
}

void print_results(size_t const fp_count, config const & cfg)
{
    double const fpr = (cfg.number_of_kmers > cfg.elements)
                         ? static_cast<double>(fp_count) / (cfg.number_of_kmers - cfg.elements)
                         : 0.0;
    std::cout << "fp_count: " << fp_count << '\n' //
              << "fp_rate: " << std::fixed << std::setprecision(3) << fpr << '\n';
}

void single_tb(config const & cfg)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{1u},
                                              seqan::hibf::bin_size{seqan::hibf::build::bin_size_in_bits(
                                                  {.fpr = cfg.fpr, .hash_count = cfg.hash, .elements = cfg.elements})},
                                              seqan::hibf::hash_function_count{cfg.hash}};
    auto agent = ibf.membership_agent();

    // Generate elements many random kmer values.
    robin_hood::unordered_set<uint64_t> inserted_values{};
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> distrib(0ULL, cfg.number_of_kmers - 1u);

    for (; inserted_values.size() < cfg.elements;)
        inserted_values.emplace(distrib(gen));

    for (uint64_t const value : inserted_values)
        ibf.emplace(value, seqan::hibf::bin_index{0u});

    // Check all possible kmer values.
    size_t fp_count{};
    for (uint64_t value{}; value < cfg.number_of_kmers; ++value)
    {
        auto & result = agent.bulk_contains(value);
        fp_count += result[0];
    }

    // There are cfg.elements many TP.
    fp_count -= cfg.elements;

    print_results(fp_count, cfg);
}

void multiple_tb(config const & cfg, size_t const bin_size)
{
    seqan::hibf::interleaved_bloom_filter ibf{seqan::hibf::bin_count{cfg.splits},
                                              seqan::hibf::bin_size{bin_size},
                                              seqan::hibf::hash_function_count{cfg.hash}};
    auto agent = ibf.membership_agent();

    // Generate elements many random kmer values.
    robin_hood::unordered_set<uint64_t> all_values{};
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> distrib(0ULL, cfg.number_of_kmers - 1u);

    for (; all_values.size() < cfg.elements;)
        all_values.emplace(distrib(gen));

    // Distribute across all bins.
    size_t counter{};
    for (uint64_t const value : all_values)
        ibf.emplace(value, seqan::hibf::bin_index{counter++ / cfg.split_elements_per_bin});

    // Check all possible kmer values.
    size_t fp_count{};
    for (uint64_t value{}; value < cfg.number_of_kmers; ++value)
    {
        auto & result = agent.bulk_contains(value);
        // A hit in any of the split bins is a FP. We only count one.
        size_t local_fp{};
        for (size_t i{}; i < cfg.splits; ++i)
        {
            local_fp += result[i];
        }
        fp_count += local_fp > 0u;
    }

    // There are cfg.elements many TP.
    fp_count -= cfg.elements;

    print_results(fp_count, cfg);
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"fpr_quality", argc, argv, sharg::update_notifications::off};
    parser.info.author = "Enrico Seiler";
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.short_description = "Inserts a given amount of k-mers into an IBF and queries all possible k-mers. "
                                    "Reports the resulting FPR for both single and split bins.";
    config cfg{};
    init_parser(parser, cfg);
    parser.parse();

    cfg.number_of_kmers = (1ULL << (2 * cfg.kmer_size));

    if (cfg.elements > cfg.number_of_kmers)
    {
        std::cout << "[WARNING] Inserting more elements than there are possible k-mers. "
                  << "Setting number of elements to number of possible k-mers.\n\n";
        cfg.elements = cfg.number_of_kmers;
    }

    cfg.split_elements_per_bin = (cfg.elements + cfg.splits - 1) / cfg.splits; // ceil for positive integers

    std::cout << "kmer: " << cfg.kmer_size << '\n';
    std::cout << "elements: " << cfg.elements << '\n';
    std::cout << "splits: " << cfg.splits << '\n';
    std::cout << "hash: " << cfg.hash << '\n';
    std::cout << "fpr: " << cfg.fpr << "\n\n";

    std::cout << "=== Single bin ===\n";
    single_tb(cfg);

    std::cout << "=== Split into " << cfg.splits << " bins ===\n";
    multiple_tb(cfg, split_bin_size_in_bits(cfg));

    std::cout << "=== Split into " << cfg.splits << " bins corrected ===\n";
    double const fpr_correction = seqan::hibf::layout::compute_fpr_correction(
        {.fpr = cfg.fpr, .hash_count = cfg.hash, .t_max = cfg.splits})[cfg.splits];
    multiple_tb(cfg, std::ceil(split_bin_size_in_bits(cfg) * fpr_correction));
}
