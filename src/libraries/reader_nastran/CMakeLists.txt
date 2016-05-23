## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(reader_nastran ${CMAKE_CURRENT_SOURCE_DIR}/src/reader_nastran.cpp)
target_include_directories(reader_nastran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(reader_nastran Triangle)
target_link_libraries(reader_nastran TriangleContainer)
target_link_libraries(reader_nastran LabelContainer)


