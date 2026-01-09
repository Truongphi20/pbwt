cmake_minimum_required(VERSION 3.15)

project(pbwt
    VERSION 2.1.0
    LANGUAGES C CXX
)

# Use a modern C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Get git commit hash (if .git exists)
set(PBWT_COMMIT_HASH "")

if (EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(
        COMMAND git describe --always --long --dirty
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        OUTPUT_VARIABLE PBWT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
endif()

# Generate version.h
configure_file(
    ${PROJECT_SOURCE_DIR}/include/version.h.in
    ${PROJECT_BINARY_DIR}/generated/version.h
    @ONLY
)

# PBWT library
add_library(pbwt_lib 
    src/pbwtMain.c
    src/array.c
    src/dict.c
    src/hash.c
    src/pbwtCore.c
    src/pbwtGeneticMap.c
    src/pbwtHtslib.c
    src/pbwtIO.c
    src/pbwtImpute.c
    src/pbwtLikelihood.c
    src/pbwtMatch.c
    src/pbwtMerge.c
    src/pbwtPaint.c
    src/pbwtSample.c
    src/utils.c
)
target_include_directories(pbwt_lib
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/generated
)

option(TEST "DO GTEST" OFF)
if (TEST)
  add_subdirectory(test)
else()
    # Build executable
    add_executable(pbwt src/main.c)

    # Linking
    target_link_libraries(pbwt 
            pbwt_lib 
            pthread
            hts
            z
            m
            bz2
            lzma
            curl
    )

    # Compiler warnings
    target_compile_options(pbwt
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
    )
endif()