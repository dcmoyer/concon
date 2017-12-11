project(c3 CXX)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-gdwarf-2 -pg -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -fopenmp")

#otherwise -rdynamic is set, linking EVERYTHING from the original libraries.
#SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
#consider -ffast-math
#add_definitions(-std=c++11)

if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "-std=gnu++11 ${CMAKE_CXX_FLAGS}")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

## setup output directories
if (CUSTOM_OUTPUT_PREFIX)
  message(STATUS "custom output prefix set to " ${CUSTOM_OUTPUT_PREFIX})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CUSTOM_OUTPUT_PREFIX}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CUSTOM_OUTPUT_PREFIX}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CUSTOM_OUTPUT_PREFIX}/bin)
else()
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

# additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#if(CUSTOM_PREFIX)
#  message(STATUS "custom prefix set to " ${CUSTOM_PREFIX})
#  set(CMAKE_LIBRARY_PATH "${CUSTOM_PREFIX}/lib")
#  set(CMAKE_LIBRARY_PATH "${CUSTOM_PREFIX}/lib64")
#  set(CMAKE_INCLUDE_PATH "${CUSTOM_PREFIX}/include")
#endif()

#
# third party
#

find_package(OpenMP)
if (NOT OPENMP_FOUND)
    message(ERROR " OpenMP not found")
endif()

#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# setup Armadillo

#find_package(Armadillo REQUIRED)
#find_library(ARMADILLO_LIB NAMES libarmadillo.a armadillo PATHS "${CMAKE_LIBRARY_PATH}")
#message(STATUS "armadillo linked from " ${ARMADILLO_LIB})
#
#find_path(ARMADILLO_INC_DIR NAMES armadillo PATHS "${CMAKE_INCLUDE_PATH}")
#include_directories(${ARMADILLO_INC_DIR})
#message(STATUS "armadillo include path " ${ARMADILLO_INC_DIR})

#link_directories(${CMAKE_LIBRARY_PATH})
#link_directories(/usr/local/fftw)

#
# local third party
#

file(GLOB sh_sources third_party/spherical_harmonics/sh/spherical_harmonics.cc)
file(GLOB sh_headers third_party/spherical_harmonics/sh/*.h)
add_library(spherical_harmonics ${sh_sources} ${sh_headers})

#
# collect sources
#

file(GLOB MeshLib_sources src/MeshLib/*.cpp)
file(GLOB MeshLib_headers src/MeshLib/*.hpp src/ops/*.h)

file(GLOB utils_sources src/utils/*.cpp)
file(GLOB utils_headers src/utils/*.hpp src/ops/*.h)

file(GLOB trk_sources src/trk/*.cpp)
file(GLOB trk_headers src/trk/*.hpp src/ops/*.h)

file(GLOB c3_sources src/c3/*.cpp src/c3/s_ops/*.cpp src/c3/g_ops/*.cpp)
file(GLOB c3_headers src/c3/*.hpp src/c3/s_ops/*.hpp src/c3/g_ops/*.hpp)
list(REMOVE_ITEM c3_headers ${CMAKE_CURRENT_SOURCE_DIR}/src/c3/prototype_sop.hpp)
list(REMOVE_ITEM c3_sources ${CMAKE_CURRENT_SOURCE_DIR}/src/c3/prototype_sop.cpp)
list(REMOVE_ITEM c3_sources src/c3/main.cpp)

file(GLOB sources src/*.cpp src/struct/*.cpp src/dist_funcs/*.cpp src/ops/*.cpp)
file(GLOB headers src/*.hpp src/struct/*.hpp src/dist_funcs/*.hpp src/ops/*.hpp)

#
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

# compile and link
add_library(MeshLib ${MeshLib_sources} ${MeshLib_headers})
add_library(utils ${utils_sources} ${utils_headers})
add_library(trk ${trk_sources} ${trk_headers})
add_library(c3 ${c3_sources} ${c3_headers})
#add_executable(viewer ${viewer_sources} ${viewer_headers})

add_executable(c3_main src/c3/main.cpp)
target_link_libraries(c3_main c3 MeshLib utils trk spherical_harmonics)
#message(STATUS "custom output prefix set to " ${CUSTOM_OUTPUT_PREFIX})
#target_link_libraries(gspn ${OPEN_BLAS} ${ARMADILLO_LIB})
#target_link_libraries(gspn ${ARMADILLO_LIBRARIES})


#target_link_libraries(MeshLib ${OPENMESH_LIBRARIES})
#target_link_libraries(viewer ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${OPENMESH_LIBRARIES})
#target_link_libraries(sd ${OPENMESH_LIBRARIES} MeshLib ${ARMADILLO_LIBRARIES})

#target_compile_features(MeshLib PRIVATE cxx_range_for)
#target_compile_features(sd PRIVATE cxx_range_for)

option(BUILD_TESTS "Build test executables" OFF)

if(BUILD_TESTS)
  add_subdirectory(test)
endif()








