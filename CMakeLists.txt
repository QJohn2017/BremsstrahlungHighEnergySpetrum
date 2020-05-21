project(BREMSSPECTRUM CXX)
cmake_minimum_required(VERSION 3.6)

set( CMAKE_VERBOSE_MAKEFILE off )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++11" )
set( BREMSSPECTRUM_INCLUDES )
set( BREMSSPECTRUM_LIBRARIES )

set (CMAKE_CXX_STANDARD 11)
cmake_policy(SET CMP0048 NEW)

# include
include_directories( "${CMAKE_SOURCE_DIR}" )

# GSL (external)
find_package(GSL)
if(GSL_FOUND)
    set(HAVE_LIBGSL 1)
    list(APPEND BREMSSPECTRUM_INCLUDES ${GSL_INCLUDE_DIR})
    list(APPEND BREMSSPECTRUM_LIBRARIES ${GSL_LIBRARIES})
else(GSL_FOUND)
    message(FATAL_ERROR "Require GSL: set GSL_HOME")
endif(GSL_FOUND)

# add source files
set( BREMSSPECTRUM_SRCS
	src/bremsstrahlung_sigma.cpp
	src/bremsstrahlung_functions.cpp
	)

# make library
include_directories (include ${BREMSSPECTRUM_INCLUDES})
add_library( BREMSSPECTRUM_LIB ${BREMSSPECTRUM_SRCS} ${BREMSSPECTRUM_INCLUDES} )

# make executables
add_executable(spectrum src/main.cpp)
target_link_libraries (spectrum BREMSSPECTRUM_LIB ${BREMSSPECTRUM_LIBRARIES})

# show summary
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_C_COMPILER ...... = ${CMAKE_C_COMPILER}")
message (STATUS "CMAKE_C_FLAGS ......... = ${CMAKE_C_FLAGS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "BREMSSPECTRUM_INCLUDES ...... = ${BREMSSPECTRUM_INCLUDES}")
message (STATUS "BREMSSPECTRUM_LIBRARIES ..... = ${BREMSSPECTRUM_LIBRARIES}")
