## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(Point3DContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/Point3DContainer.cpp)
target_include_directories(Point3DContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Point3D ${CMAKE_CURRENT_SOURCE_DIR}/src/Point3D.cpp)
target_include_directories(Point3D PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)



