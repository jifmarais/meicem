## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

# Math libraries
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})


add_library(Vector ${CMAKE_CURRENT_SOURCE_DIR}/src/Vector.cpp)
target_include_directories(Vector PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(ComplexMatrix ${CMAKE_CURRENT_SOURCE_DIR}/src/ComplexMatrix.cpp)
target_include_directories(ComplexMatrix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ComplexMatrix ${ARMADILLO_LIBRARIES})

add_library(Quadrature ${CMAKE_CURRENT_SOURCE_DIR}/src/Quadrature.cpp)
target_include_directories(Quadrature PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Quadrature ${ARMADILLO_LIBRARIES})
target_link_libraries(Quadrature Node)
target_link_libraries(Quadrature Triangle)
