## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(Vector ${CMAKE_CURRENT_SOURCE_DIR}/src/Vector.cpp)
target_include_directories(Vector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)



