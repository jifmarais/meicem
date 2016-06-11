## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/math)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/math/test)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mesh)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mesh/test)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/reader_nastran)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/reader_nastran/test)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/reader_wilco_input)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/reader_wilco_input/test)



