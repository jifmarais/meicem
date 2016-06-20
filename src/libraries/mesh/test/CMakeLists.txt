## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_executable(mesh.test ${CMAKE_CURRENT_SOURCE_DIR}/src/main.test.cpp)

# BOOST libraries
find_package( Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries(mesh.test ${Boost_LIBRARIES})
target_link_libraries(mesh.test NodeContainer)
target_link_libraries(mesh.test Node)
target_link_libraries(mesh.test Triangle)
target_link_libraries(mesh.test TriangleContainer)
target_link_libraries(mesh.test LabelContainer)
target_link_libraries(mesh.test Edge)
target_link_libraries(mesh.test EdgeContainer)

target_include_directories(mesh.test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)



