## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_executable(reader_wilco_input.test ${CMAKE_CURRENT_SOURCE_DIR}/src/main.test.cpp)

# BOOST libraries
find_package( Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries(reader_wilco_input.test ${Boost_LIBRARIES})
target_link_libraries(reader_wilco_input.test reader_wilco_input)
target_link_libraries(reader_wilco_input.test NodeContainer)
target_link_libraries(reader_wilco_input.test Node)
#target_link_libraries(reader_wilco_input.test Triangle)
target_link_libraries(reader_wilco_input.test TriangleContainer)

#target_include_directories(reader_nastran.test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)



