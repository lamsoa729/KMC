cmake_minimum_required(VERSION 3.7)
project(KMC)

# Set binary output directory to build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(ALLOW_DUPLICATE_CUSTOM_TARGETS TRUE)

if(DEBUG OR UNIT_TESTS)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

if(CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
    message("Building KMC in Debug mode")
else()
  # otherwise set build to release mode by default
    set(CMAKE_BUILD_TYPE "Release")
    message("Building KMC in Release mode")
endif()

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory \
with a CMakeLists.txt file). Please make a build subdirectory. Feel free to \
remove CMakeCache.txt and CMakeFiles.")
endif()

# Add executables
add_subdirectory("KMC")
