cmake_minimum_required(VERSION 3.14)
project(KATCH)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_VERBOSE_MAKEFILE ON)

#set(CMAKE_CXX_FLAGS "-DUSE_PARALLELISM")
set(CMAKE_CXX_FLAGS_DEBUG  "-g -ggdb -O0 -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/../bin)

if(UNIX AND NOT APPLE)
    message(STATUS ">>> Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include -fopenmp")
endif ()

if(APPLE)
    message(STATUS ">>> Apple")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")
    include_directories("/usr/local/Cellar/libomp/10.0.0/include")
    include_directories("/usr/local/Cellar/google-sparsehash/2.0.4/include")
    include_directories("/usr/local/Cellar/boost/1.75.0_1/include")
endif()


include_directories(./)
include_directories(osm_analyzer)
set(SRC

        datastr/base/double.h
        datastr/base/interval.h
        datastr/base/point.h
        datastr/base/pwl_ttf.h
        datastr/base/ttf_wrapper.h
        datastr/base/bit_vector.cpp
        datastr/base/bit_vector.h
        datastr/base/emulate_gcc_builtin.h
        datastr/base/spline.h

        datastr/graph/basic.h
        datastr/graph/dynamic_graph.h
        datastr/graph/order_edge_data.h
        datastr/graph/search_graph.h
        datastr/graph/dynamic_search_graph.h
        datastr/graph/static_search_graph.h
        datastr/graph/forward_search_graph.h
        datastr/graph/static_forward_search_graph.h


        datastr/cpd/cpd.h
        datastr/cpd/cpd.cpp
        datastr/cpd/range.h
        datastr/cpd/binary_search.h
        datastr/cpd/bw_cpd.cpp
        datastr/cpd/bw_cpd.h
        datastr/cpd/rev_cpd.h
        datastr/cpd/rev_bw_cpd.h
        datastr/cpd/rev_bw_cpd.cpp

        io/binary.h
        io/btch_format.h
        io/demands_format.h
        io/edge_info.h
        io/tpgr_format.h
        io/vec_io.h
        io/tpgr_format_noise.h

        prepro/edge_description.h
        prepro/ordering.h
        prepro/witness_search.h


        query/BTCH.h
        query/BSTCH.h
        query/FTCH_TCPD.h
        query/STCH_TCPD.h
        query/FTCH_L.h
        query/BMTCH.h
        query/BMTCH_query.h
        query/FTCH_CPD.h
        query/FMTCH_TCPD.h
        query/FMTCH_TCPD_query.h
        query/FTCH_RTPD.h
        query/BTCH_L.h

        cpd_prepro/tch_downward_dijkstra_query.h
        cpd_prepro/st_dijkstra_query.h
        cpd_prepro/tch_dijkstra_query.h
        cpd_prepro/tch_updown_query.h
        cpd_prepro/reverse_tch_updown_query.h
        cpd_prepro/tch_upward_dijkstra_query.h

        util/id_queue.h
        util/misc.h
        util/my_timer.h
        util/my_timer.cpp

        context/search_context.h
        context/bi_search_context.h
        )

add_executable(run_ordering ${SRC} experiments/run_ordering.cpp)
add_executable(run_experiments ${SRC} experiments/run_experiments.cpp)
add_executable(generate_queries ${SRC} experiments/generate_queries.cpp)
add_executable(split_data ${SRC} experiments/split_data.cpp)
add_executable(construct_fw_tch_cpd ${SRC} experiments/construct_fw_tch_cpd.cpp)
add_executable(construct_bw_tch_cpd ${SRC} experiments/construct_bw_tch_cpd.cpp)
add_executable(construct_stch_mtch ${SRC} experiments/construct_stch_mtch.cpp)
add_executable(generate_landmarks ${SRC} experiments/generate_landmarks.cpp )
add_executable(construct_reverse_fw_tch_cpd ${SRC} experiments/construct_reverse_fw_tch_cpd.cpp)

find_library(omp_lib omp "/usr/local/opt")

if (omp_lib)
    target_link_libraries(run_ordering ${omp_lib})
    target_link_libraries(construct_stch_mtch ${omp_lib})
    target_link_libraries(construct_fw_tch_cpd ${omp_lib})
    target_link_libraries(construct_bw_tch_cpd ${omp_lib})
    target_link_libraries(construct_reverse_fw_tch_cpd ${omp_lib})
endif()
