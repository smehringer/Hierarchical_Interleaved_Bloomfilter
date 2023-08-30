#include <gtest/gtest.h> // for Test, TestInfo, Message, TestPartResult, TEST, EXPECT_EQ

#include <cstddef>   // for size_t
#include <stdexcept> // for runtime_error
#include <vector>    // for allocator, vector

#include <hibf/layout/data_store.hpp>     // for data_store
#include <hibf/layout/layout.hpp>         // for layout
#include <hibf/layout/simple_binning.hpp> // for simple_binning
#include <hibf/test/expect_range_eq.hpp>  // for expect_range_eq, EXPECT_RANGE_EQ

TEST(simple_binning_test, small_example)
{
    seqan::hibf::layout::layout hibf_layout;
    std::vector<size_t> kmer_counts{100, 40, 20, 20};

    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout,
                                 .kmer_counts = &kmer_counts,
                                 .fpr_correction = std::vector<double>(65, 1.0)};

    seqan::hibf::layout::simple_binning algo{data, 9};
    size_t max_bin = algo.execute();

    std::vector<seqan::hibf::layout::layout::user_bin> expected{{3, {}, 1, 0},
                                                                {2, {}, 1, 1},
                                                                {1, {}, 2, 2},
                                                                {0, {}, 5, 4}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, uniform_distribution)
{
    seqan::hibf::layout::layout hibf_layout;
    std::vector<size_t> kmer_counts{20, 20, 20, 20};

    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout,
                                 .kmer_counts = &kmer_counts,
                                 .fpr_correction = std::vector<double>(65, 1.0)};

    seqan::hibf::layout::simple_binning algo{data, 4u};
    size_t max_bin = algo.execute();

    std::vector<seqan::hibf::layout::layout::user_bin> expected{{3, {}, 1, 0},
                                                                {2, {}, 1, 1},
                                                                {1, {}, 1, 2},
                                                                {0, {}, 1, 3}};

    EXPECT_RANGE_EQ(hibf_layout.user_bins, expected);
    EXPECT_EQ(max_bin, 0u);
}

TEST(simple_binning_test, user_bins_must_be_smaller_than_technical_bins)
{
    seqan::hibf::layout::layout hibf_layout;

    std::vector<size_t> kmer_counts{100, 40, 20, 20};

    seqan::hibf::layout::data_store data{.hibf_layout = &hibf_layout,
                                 .kmer_counts = &kmer_counts,
                                 .fpr_correction = std::vector<double>(65, 1.0)};

    EXPECT_THROW((seqan::hibf::layout::simple_binning{data, 2}), std::runtime_error);
}
