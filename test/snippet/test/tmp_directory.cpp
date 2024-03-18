// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h> // for Test, AssertionResult, TestInfo, EXPECT_TRUE, Message, TEST, TestPar...

#include <filesystem> // for remove
#include <fstream>    // for char_traits, operator<<, basic_ostream, ofstream

#include <hibf/test/sandboxed_path.hpp> // for operator/, sandboxed_path
#include <hibf/test/tmp_directory.hpp>  // for tmp_directory

TEST(snippet_tmp_directory, tmp_directory_)
{
    // create a directory folder
    seqan::hibf::test::tmp_directory tmp{};

    // Some function that should creates temporary files and removes them again
    {
        std::ofstream ofs{tmp.path() / "somefile.txt"};
        ofs << "Hello World!";
        ofs.close();

        std::filesystem::remove(tmp.path() / "somefile.txt");
    }

    // check that everything was cleanup properly
    EXPECT_TRUE(tmp.empty());
}
