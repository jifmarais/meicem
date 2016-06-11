## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_executable(math.test ${CMAKE_CURRENT_SOURCE_DIR}/src/main.test.cpp)

# BOOST libraries
find_package( Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries(math.test ${Boost_LIBRARIES})
target_link_libraries(math.test Vector)

target_include_directories(math.test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)



