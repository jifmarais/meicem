## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_executable(mom.test ${CMAKE_CURRENT_SOURCE_DIR}/src/main.test.cpp)

# BOOST libraries
find_package( Boost REQUIRED COMPONENTS unit_test_framework)
target_link_libraries(mom.test ${Boost_LIBRARIES})
target_link_libraries(mom.test MoM)
target_link_libraries(mom.test NearFieldValue)
target_link_libraries(mom.test NearFieldContainer)
target_link_libraries(mom.test PlaneWave)

target_include_directories(mom.test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)



