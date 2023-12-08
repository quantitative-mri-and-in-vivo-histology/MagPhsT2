# CMake compatibility issues
CMAKE_MINIMUM_REQUIRED( VERSION 3.0 )
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# Project name
PROJECT(EPG)

# Generate compile_commands.json for autocompletion
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Installation prefix
SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Compiler flags
SET(CMAKE_CXX_FLAGS "-Wall -O3 -g -fPIC -std=c++11")
MESSAGE("-- CXX Flags: " ${CMAKE_CXX_FLAGS})

# Needed for pyepg.cxx in build directory during make
INCLUDE_DIRECTORIES( ${EPG_SOURCE_DIR})
MESSAGE("-- EPG include dir: " ${EPG_SOURCE_DIR})

# Add CNPY to read and write python arrays in npy format
ADD_SUBDIRECTORY(external/cnpy)

# Libs to create
ADD_LIBRARY(epgtools epg.cpp epg.h)
SET_TARGET_PROPERTIES(epgtools PROPERTIES OUTPUT_NAME epg)

# Name of executable
ADD_EXECUTABLE(epg main.cpp)

# # Link it
TARGET_LINK_LIBRARIES(epg epgtools cnpy)

MESSAGE("-- Installation directory: " ${CMAKE_INSTALL_PREFIX})

# Install header files
INSTALL(FILES epg.h DESTINATION include)

# Install additional targets (binaries, libs)
INSTALL(TARGETS epg epgtools cnpy
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
