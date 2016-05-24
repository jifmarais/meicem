## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(NodeContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/NodeContainer.cpp)
target_include_directories(NodeContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Node ${CMAKE_CURRENT_SOURCE_DIR}/src/Node.cpp)
target_include_directories(Node PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(TriangleContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/TriangleContainer.cpp)
target_include_directories(TriangleContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Triangle ${CMAKE_CURRENT_SOURCE_DIR}/src/Triangle.cpp)
target_include_directories(Triangle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(LabelContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/LabelContainer.cpp)
target_include_directories(LabelContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)



