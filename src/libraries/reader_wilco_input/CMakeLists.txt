## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(reader_wilco_input ${CMAKE_CURRENT_SOURCE_DIR}/src/reader_wilco_input.cpp)
target_include_directories(reader_wilco_input PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(reader_wilco_input Triangle)
target_link_libraries(reader_wilco_input TriangleContainer)
target_link_libraries(reader_wilco_input LabelContainer)


