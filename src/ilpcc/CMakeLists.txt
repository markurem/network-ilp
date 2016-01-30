cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(ilpcc)
include_directories(include)
file(GLOB headers_insegt
    include/andres/*
    include/andres/ilp/*
)

##############################################################################
# MSVC-specific settings
##############################################################################
if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

##############################################################################
# OpenMP
##############################################################################
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found")
else()
    message("OpenMP not found")
endif()

##############################################################################
# packages
##############################################################################
find_package(CPLEX)
if(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    include_directories(${CONCERT_INCLUDE_DIR})
    add_definitions(-DILOUSESTL)
    add_definitions(-DIL_STD)
endif()

##############################################################################
# targets
##############################################################################
if(CPLEX_FOUND)
    add_executable(test-ilpcc src/unit-tests/ilpcc.cxx ${headers_insegt})
    target_link_libraries(test-ilpcc
        ${CMAKE_THREAD_LIBS_INIT}
        ${ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY} ${CONCERT_LIBRARY}
    )
endif()

