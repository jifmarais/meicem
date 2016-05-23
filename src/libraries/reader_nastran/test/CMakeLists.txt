## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_executable(reader_nastran.test ${CMAKE_CURRENT_SOURCE_DIR}/src/main.test.cpp)

# BOOST libraries
find_package( Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries(reader_nastran.test ${Boost_LIBRARIES})
target_link_libraries(reader_nastran.test reader_nastran)
target_link_libraries(reader_nastran.test Point3DContainer)
target_link_libraries(reader_nastran.test Point3D)
#target_link_libraries(reader_nastran.test Triangle)
target_link_libraries(reader_nastran.test TriangleContainer)

#target_include_directories(reader_nastran.test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)



