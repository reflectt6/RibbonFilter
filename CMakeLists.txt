cmake_minimum_required(VERSION 3.22.1)
project(ribbon_filter)

set(CMAKE_CXX_STANDARD 17)

include_directories(./
        $ENV{JAVA_HOME}/include/
        $ENV{JAVA_HOME}/include/darwin)

#set(SOURCE
#        util/hash.cc
#        util/murmurhash.cc
#        util/ribbon_config.cc
#        util/xxhash.cc
#        library.cpp ribbon-invoke/ribbon_compare_bloom_prune.cc)

set(TESTS
        ribbon-lib/hash.cc
        ribbon-lib/xxhash.cc
        ribbon-lib/ribbon_config.cc
#        ribbon-invoke/org_example_JNIAdvancedDemo.cpp
#        ribbon-invoke/ribbon_native_invoke_simple_demo.cc
        ribbon-invoke/org_apache_hadoop_hbase_util_RibbonHelper.cpp
        )

add_library(ribbon_filter SHARED ${TESTS})

#add_executable(ribbon_test ${TESTS})
