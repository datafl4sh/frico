cmake_minimum_required(VERSION 3.20)
project(mommy)

set(CMAKE_CXX_STANDARD 23)

set (CMAKE_MODULE_PATH
     ${CMAKE_MODULE_PATH}
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set (CMAKE_CXX_FLAGS "-fsanitize=address")
#set (CMAKE_CXX_FLAGS "-g")

find_package(GMSH REQUIRED)
set(LINK_LIBS ${LINK_LIBS} GMSH)

find_package(SILO REQUIRED)
set(LINK_LIBS ${LINK_LIBS} SILO)

find_package(Eigen3 REQUIRED)
set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

if (NOT DEFINED HIGHFIVE_USE_EIGEN)
    set(HIGHFIVE_USE_EIGEN On)
endif()

if (NOT DEFINED HIGHFIVE_USE_BOOST)
    set(HIGHFIVE_USE_BOOST Off)
endif()

add_subdirectory(deps/HighFive)
set(LINK_LIBS ${LINK_LIBS} HighFive)

set(SOURCES main.cpp geom_mesh.cpp input_gmsh.cpp output_silo.cpp quadratures.cpp)

add_executable(mommy ${SOURCES})

target_link_libraries(mommy PUBLIC
    ${LINK_LIBS}
)
