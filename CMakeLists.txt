## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

# Compiler options
if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++11 ")
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

# Local libraries
target_link_libraries(meicem Point3DContainer)
target_link_libraries(meicem Point3D)
#target_link_libraries(meicem PointContainer)
#target_link_libraries(meicem Triangle)


# Tests
enable_testing ()
add_test (NAME mesh COMMAND mesh.test)
add_test (NAME reader_nastran COMMAND reader_nastran.test)

