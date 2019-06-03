cmake_minimum_required(VERSION 3.12.4)
project(kmcx)

# Set binary output directory to build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

#set(GXX_COVERAGE_COMPILE_FLAGS "-Wunused-value -Wreturn-type")
#add_definitions(${GXX_COVERAGE_COMPILE_FLAGS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/.CMake_FFTW")

if(DEBUG)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

if(TESTS)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

if(CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
  # if debug mode, do tests
  if(TESTS)
    message("Building unit tests")
  else()
    message("Building kmcx in Debug mode")
  endif()
else()
  # otherwise set build to release mode by default
  set(CMAKE_BUILD_TYPE "Release")
  message("Building simcore in Release mode")
endif()

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Add executables

# Reset build output directory to default

# Configure unit tests
if (TESTS)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  enable_testing()
  add_subdirectory("tests")
else()
  add_subdirectory("src")
endif(TESTS)

# Build the documentation
# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  # set input and output files
  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

  # configure the file
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message("Adding Doxygen build target")

  # Can build docs-only with the `make docs` target
  add_custom_target( docs
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating documentation with Doxygen"
      VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the documentation")
endif (DOXYGEN_FOUND)
