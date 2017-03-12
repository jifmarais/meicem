## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

# Compiler options
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -std=c++11 -O2 ")
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++11 ")
endif(UNIX)

# Get some debug info
SET(CMAKE_VERBOSE_MAKEFILE on )
set(CTEST_OUTPUT_ON_FAILURE TRUE )

#CONFIGURATION OF THE PROJECT
#NAME OF THE PROJECT
project(meicem)

#OUTPUT OF THE PROJECT
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE})

# Add the executable
add_executable(meicem ${CMAKE_CURRENT_SOURCE_DIR}/src/meicem.cpp )
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/libraries)

# BOOST libraries
find_package( Boost REQUIRED COMPONENTS program_options )
target_link_libraries(meicem ${Boost_LIBRARIES})

# Math libraries
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(meicem ${ARMADILLO_LIBRARIES})

## Math libraries - LAPACK
#find_package( LAPACK REQUIRED )
#include_directories(${LAPACK_INCLUDE_DIR})
#target_link_libraries(meicem ${LAPACK_LIBRARIES})

## Math libraries - BLAS
#find_package( BLAS REQUIRED )
#include_directories(${BLAS_INCLUDE_DIR})
#target_link_libraries(meicem ${BLAS_LIBRARIES})

# Local libraries
#target_link_libraries(meicem Vector)
target_link_libraries(meicem Quadrature)
target_link_libraries(meicem EdgeContainer)
target_link_libraries(meicem NodeContainer)
target_link_libraries(meicem NearFieldContainer)
target_link_libraries(meicem NearFieldValue)
target_link_libraries(meicem TriangleContainer)
target_link_libraries(meicem ComplexMatrix)
target_link_libraries(meicem reader_wilco_input)
target_link_libraries(meicem reader_nastran)
target_link_libraries(meicem MoM)
#target_link_libraries(meicem Node)
#target_link_libraries(meicem Triangle)


# Tests
enable_testing ()
add_test (NAME math COMMAND math.test)
add_test (NAME mesh COMMAND mesh.test)
add_test (NAME reader_nastran COMMAND reader_nastran.test)
add_test (NAME reader_wilco_reader COMMAND reader_wilco_input.test)
add_test (NAME mom_3d COMMAND mom.test)

