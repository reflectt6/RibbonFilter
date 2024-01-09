cmake_minimum_required(VERSION 3.24)
project(ribbon_filter)

set(CMAKE_CXX_STANDARD 17)

include_directories(ribbon-lib)

#set(SOURCE
#        util/hash.cc
#        util/murmurhash.cc
#        util/ribbon_config.cc
#        util/xxhash.cc
#        library.cpp ribbon-invoke/ribbon_compare_bloom_prune.cc)

set(TESTS
        ribbon-lib/util/hash.cc
        ribbon-lib/util/murmurhash.cc
        ribbon-lib/util/xxhash.cc
        ribbon-lib/util/ribbon_config.cc
        ribbon-invoke/ribbon_native_invoke_demo.cc
        )

#add_library(ribbon_filter SHARED ${SOURCE})

add_executable(ribbon_test ${TESTS})
