## Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.5)

add_library(MoM ${CMAKE_CURRENT_SOURCE_DIR}/src/MoM.cpp)
target_include_directories(MoM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(MoM TriangleContainer)
target_link_libraries(MoM EdgeContainer)
target_link_libraries(MoM ComplexMatrix)
target_link_libraries(MoM Quadrature)
target_link_libraries(MoM PlaneWave)
target_link_libraries(MoM EMconst)

add_library(NearFieldContainer ${CMAKE_CURRENT_SOURCE_DIR}/src/NearFieldContainer.cpp)
target_include_directories(NearFieldContainer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(NearFieldContainer NodeContainer)
target_link_libraries(NearFieldContainer ComplexMatrix)

add_library(NearFieldValue ${CMAKE_CURRENT_SOURCE_DIR}/src/NearFieldValue.cpp)
target_include_directories(NearFieldValue PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(PlaneWave ${CMAKE_CURRENT_SOURCE_DIR}/src/PlaneWave.cpp)
target_include_directories(PlaneWave PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(PlaneWave Vector)
target_link_libraries(PlaneWave EMconst)
target_link_libraries(PlaneWave Node)

