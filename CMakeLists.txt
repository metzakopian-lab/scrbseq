cmake_minimum_required(VERSION 3.10)
project(SplitFiles)


find_package(Boost 1.66.0 COMPONENTS iostreams system)
add_executable(splitFQ split_fq.cpp)
add_executable(splitSAM split_sam.cpp)


target_link_libraries(splitFQ ${Boost_LIBRARIES})

target_link_libraries(splitSAM ${Boost_LIBRARIES})
