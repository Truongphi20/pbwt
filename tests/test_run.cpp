#include <gtest/gtest.h>
#include "pbwtMain.h"

#include <vector>
#include <string>
#include <cstring>

#ifndef TEST_DATA_FOLDER
#define TEST_DATA_FOLDER
#endif

TEST(DRY_RUN, test_simple_data){
    std::vector<std::string> args{
        "pbwt",
        "-checkpoint", "10000",
        "-readVcfGT",  std::string(TEST_DATA_FOLDER) + "/OMNI.vcf",
        "-writeAll", "data"
    };

    std::vector<char*> argv;
    argv.reserve(args.size());

    for (const auto& s : args) {
        char* buf = new char[s.size() + 1];
        std::strcpy(buf, s.c_str());
        argv.push_back(buf);
    }

    int argc = static_cast<int>(argv.size());
    int ret = pbwtMain(argc, argv.data());
    
    EXPECT_EQ(ret, 0);
}