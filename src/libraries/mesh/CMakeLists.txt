## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(NodeContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/NodeContainer.cpp)
target_include_directories(NodeContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(NodeContainer ComplexMatrix)

add_library(Node ${CMAKE_CURRENT_SOURCE_DIR}/src/Node.cpp)
target_include_directories(Node PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Node ComplexMatrix)

add_library(Triangle ${CMAKE_CURRENT_SOURCE_DIR}/src/Triangle.cpp)
target_include_directories(Triangle PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Triangle ComplexMatrix)
target_link_libraries(Triangle Node)
target_link_libraries(Triangle Edge)

add_library(TriangleContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/TriangleContainer.cpp)
target_include_directories(TriangleContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(TriangleContainer ComplexMatrix)
target_link_libraries(TriangleContainer Node)
target_link_libraries(TriangleContainer Edge)

add_library(LabelContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/LabelContainer.cpp)
target_include_directories(LabelContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(Edge ${CMAKE_CURRENT_SOURCE_DIR}/src/Edge.cpp)
target_include_directories(Edge PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Edge NodeContainer)
target_link_libraries(Edge TriangleContainer)

add_library(EdgeContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/EdgeContainer.cpp)
target_include_directories(EdgeContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(EdgeContainer Edge)
target_link_libraries(EdgeContainer NodeContainer)
target_link_libraries(EdgeContainer TriangleContainer)


