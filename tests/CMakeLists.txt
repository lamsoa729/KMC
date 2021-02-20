# Force recent cmake and C++
cmake_minimum_required (VERSION 3.4)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Require Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Set libraries and include directories
set(LIB ${LIB} ${Boost_LIBRARIES})
set(INCLUDES ${INCLUDES} ${Boost_INCLUDE_DIRS})

add_executable(test_kmc test_main.cpp)

# Set the include directory for builds (tests, examples, etc) and for clients
target_include_directories(KMC INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:>)
    
target_link_libraries(test_kmc PRIVATE KMC)

add_test(NAME kmc_test COMMAND test_kmc)