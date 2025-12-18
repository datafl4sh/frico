cmake_minimum_required(VERSION 3.20)
project(mommy)

set(CMAKE_CXX_STANDARD 20)

set (CMAKE_MODULE_PATH
     ${CMAKE_MODULE_PATH}
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set (CMAKE_CXX_FLAGS "-fsanitize=address")

find_package(GMSH REQUIRED)
set(LINK_LIBS ${LINK_LIBS} GMSH)

find_package(SILO REQUIRED)
set(LINK_LIBS ${LINK_LIBS} SILO)

set(SOURCES main.cpp input_gmsh.cpp)

add_executable(mommy ${SOURCES})

target_link_libraries(mommy PUBLIC
    ${LINK_LIBS}
)
