cmake_minimum_required(VERSION 3.10)
project(SplitFiles)

#set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_STATIC_RUNTIME ON)
#set(ZLIB_USE_STATIC_LIBS ON)
#set(ZLIB_USE_STATIC_RUNTIME ON)
#find_package(ZLIB REQUIRED)
#find_package(BZIP2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS iostreams system program_options)
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


set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  
  set(CMAKE_BUILD_TYPE Release)
endif()



set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")



target_link_libraries(splitFQ ${Boost_LIBRARIES} )
target_link_libraries(splitSAM  ${Boost_LIBRARIES} )
target_link_libraries(splitFQpaired ${Boost_LIBRARIES} )



install(TARGETS splitFQ splitSAM splitFQpaired
	DESTINATION bin)
