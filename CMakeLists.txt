cmake_minimum_required(VERSION 3.20)
project(frico)

set(CMAKE_CXX_STANDARD 23)

set (CMAKE_MODULE_PATH
     ${CMAKE_MODULE_PATH}
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set (CMAKE_CXX_FLAGS "-g -fsanitize=address")
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

set(MKL_INTERFACE_FULL intel_lp64)
find_package(MKL CONFIG)
if (MKL_FOUND)
  set(LINK_LIBS ${LINK_LIBS} MKL::MKL)
  add_compile_definitions(-DHAVE_MKL)
endif()

add_subdirectory(deps/HighFive)
set(LINK_LIBS ${LINK_LIBS} HighFive)

set(SOURCES src/main.cpp src/geom_mesh.cpp src/input_gmsh.cpp
            src/output_silo.cpp src/quadratures.cpp src/rwg_basis.cpp
)

add_executable(frico ${SOURCES})

target_link_libraries(frico PUBLIC
    ${LINK_LIBS}
)
