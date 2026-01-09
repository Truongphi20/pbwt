include(FetchContent)

# Invoke Gtest
fetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(googletest)
enable_testing()

# Set up a run from the tutorial: https://choishingwan.github.io/PRS-Tutorial/
add_executable(test_run test_run.c)
target_link_libraries(test_run GTest::gtest_main 
                                       pbwt_lib 
                                        pthread
                                        hts
                                        z
                                        m
                                        bz2
                                        lzma
                                        curl
)

target_compile_options(test_run
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)

target_compile_definitions(test_run PRIVATE TEST_DATA_FOLDER="${CMAKE_SOURCE_DIR}/tests/data")
target_include_directories(test_run PUBLIC ${CMAKE_SOURCE_DIR}/include)
add_test(NAME test_run COMMAND test_run)