cmake_minimum_required(VERSION 3.12.4)
project(simcore_kmc)

# Set binary output directory to build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Build unit tests by default
option(DO_TESTS "Build unit tests" ON)

if (CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
  # if debug mode, do tests
  set(DO_TESTS ON)
else ()
  # otherwise set build to release mode by default
  set(CMAKE_BUILD_TYPE Release)
endif()

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Add executables
#add_subdirectory("src")

# Reset binary output directory to default
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
#add_subdirectory(extern/Catch2)

# Configure unit tests
if (DO_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(DO_TESTS)

# Build documentation by default
option(BUILD_DOC "Build documentation" ON)

if (CMAKE_BUILD_TYPE MATCHES "^[Dd]ebug")
    # don't build the docs in debug mode
    set(BUILD_DOC OFF)
endif()

# Build the documentation
if (BUILD_DOC)
  # check if Doxygen is installed
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      # set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      # configure the file
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("Doxygen build started")

      # Can build docs-only with the `make docs` target
      add_custom_target( docs ALL
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM )
  else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif(BUILD_DOC)
