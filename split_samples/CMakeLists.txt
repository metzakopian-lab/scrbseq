cmake_minimum_required(VERSION 3.10)
project(SplitFiles)


find_package(Boost 1.63.0 COMPONENTS iostreams system)
add_executable(splitFQ
  split_fq.cpp
  io.cpp
  FastQ.cpp
  )
add_executable(splitSAM split_sam.cpp)
add_executable(splitFQpaired
  split_fq_paired.cpp
  io.cpp
  FastQ.cpp
  )


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

target_link_libraries(splitFQ ${Boost_LIBRARIES})

target_link_libraries(splitSAM ${Boost_LIBRARIES})
target_link_libraries(splitFQpaired ${Boost_LIBRARIES})



install(TARGETS splitFQ splitSAM splitFQpaired
	DESTINATION bin)
