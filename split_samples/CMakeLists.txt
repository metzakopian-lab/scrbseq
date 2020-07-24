cmake_minimum_required(VERSION 3.10)
project(SplitFiles)


find_package(Boost 1.63.0 COMPONENTS iostreams system)
add_executable(splitFQ split_fq.cpp)
add_executable(splitSAM split_sam.cpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS OFF)


target_link_libraries(splitFQ ${Boost_LIBRARIES})

target_link_libraries(splitSAM ${Boost_LIBRARIES})



install(TARGETS splitFQ splitSAM
	DESTINATION bin)
