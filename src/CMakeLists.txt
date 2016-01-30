CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
set(VERBOSE 1)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
if(NOT CMAKE_BUILD_TYPE)
  #set(CMAKE_BUILD_TYPE "DEBUG")
  set(CMAKE_BUILD_TYPE "RELEASE")
  #set(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
  #set(CMAKE_BUILD_TYPE "MINSIZEREL")
endif()

find_package(Boost 1.54.0)
if(Boost_FOUND)
  include_directories("${Boost_INCLUDE_DIRS}" "/usr/include/python2.7")
  link_directories(${Boost_LIBRARY_DIRS})
  find_package(Boost 1.54.0 COMPONENTS python)
elseif(NOT Boost_FOUND)
  message(FATAL_ERROR "Unable to find correct Boost version. Did you set BOOST_ROOT?")
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions("-Wall")
else()
  message(FATAL_ERROR "CMakeLists.txt has not been tested/written for your compiler.")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP found")
else()
    message("OpenMP not found")
endif()

find_package(CPLEX REQUIRED)
if(CPLEX_FOUND)
    include_directories(${CPLEX_INCLUDE_DIR})
    include_directories(${CONCERT_INCLUDE_DIR})
    add_definitions(-DILOUSESTL)
    add_definitions(-DIL_STD)
endif()

include_directories(ilpcc/include)
include_directories(graph/include)

###############################################################################
## ADD TARGETS

# VESSEL NETWORK OPTIMIZER
add_library(netopt SHARED py_vnilp.cxx)
set_property(TARGET netopt PROPERTY CXX_STANDARD 11) # this requires cmake>=v3.1
target_link_libraries(netopt
        ${CMAKE_THREAD_LIBS_INIT}
        ${ILOCPLEX_LIBRARY}
        ${CPLEX_LIBRARY}
        ${CONCERT_LIBRARY}
        ${Boost_LIBRARIES}
)

###############################################################################
## python setup
add_subdirectory(python)
