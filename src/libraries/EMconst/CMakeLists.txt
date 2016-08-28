## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

#add_library(EMconst)
#add_library(EMconst ${CMAKE_CURRENT_SOURCE_DIR}/include/EMconst.hpp)
add_library(EMconst INTERFACE)
target_include_directories(EMconst INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(EMconst PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#add_custom_target(EMconst SOURCES EMconst)
