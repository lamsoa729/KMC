# Find all headers in the project include folder
file(GLOB HEADERS CONFIGURE_DEPENDS "*.hpp")

find_package(Boost REQUIRED)
set(LIB ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

# Make a static library (tradeoff of memory for better performance)
add_library(KMC INTERFACE)

# Add an alias so the target can be used while building
add_library(KMC::KMC ALIAS KMC)

# Set the include directory for builds (tests, examples, etc) and for clients
target_include_directories(KMC INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:>)

target_link_libraries(KMC ${LIB})

# Set install targets and destinations
install(TARGETS KMC
    EXPORT KMCTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/KMC)

# Export additional install files (for clients using CMake)
install(EXPORT KMCTargets
    NAMESPACE KMC::
    FILE KMCConfig.cmake
    DESTINATION lib/cmake/KMC)

