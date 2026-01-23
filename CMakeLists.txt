cmake_minimum_required(VERSION 3.20)
project(frico)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set (CMAKE_MODULE_PATH
     ${CMAKE_MODULE_PATH}
     ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set (CMAKE_CXX_FLAGS "-g -fsanitize=address")
#set (CMAKE_CXX_FLAGS "-g")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

######################################################################
## RPATH settings
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

######################################################################
## GMSH
find_package(GMSH REQUIRED)
set(LINK_LIBS ${LINK_LIBS} GMSH)

######################################################################
## SILO
find_package(SILO REQUIRED)
set(LINK_LIBS ${LINK_LIBS} SILO)

######################################################################
## Eigen
find_package(Eigen3 REQUIRED)
set(LINK_LIBS ${LINK_LIBS} Eigen3::Eigen)

######################################################################
## HighFive
if (NOT DEFINED HIGHFIVE_USE_EIGEN)
    set(HIGHFIVE_USE_EIGEN On)
endif()

if (NOT DEFINED HIGHFIVE_USE_BOOST)
    set(HIGHFIVE_USE_BOOST Off)
endif()

add_subdirectory(deps/HighFive)
set(LINK_LIBS ${LINK_LIBS} HighFive)

######################################################################
## Update submodules if necessary
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()

    execute_process( COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "DISKPP_GIT_REVISION"
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
    set(DISKPP_GIT_REVISION 0)
endif()

######################################################################
## Intel MKL: this MUST remain OPTIONAL
set(MKL_INTERFACE_FULL intel_lp64)
find_package(MKL CONFIG)
if (MKL_FOUND)
  set(LINK_LIBS ${LINK_LIBS} MKL::MKL)
  add_compile_definitions(-DHAVE_MKL)
endif()

######################################################################
## Sources and targets

option(OPT_ENABLE_OPENMP "Enable OpenMP" ON)
if (OPT_ENABLE_OPENMP)
    find_package(OpenMP)
    if (OpenMP_FOUND)
        set(LINK_LIBS ${LINK_LIBS} OpenMP::OpenMP_CXX)
    endif()
endif()


set(SOURCES src/main.cpp src/geom_mesh.cpp src/input_gmsh.cpp
            src/output_silo.cpp src/quadratures.cpp src/rwg_basis.cpp
            src/utils.cpp src/sources.cpp
)

add_executable(frico ${SOURCES})

target_link_libraries(frico PUBLIC ${LINK_LIBS})
