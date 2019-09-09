find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(LIB ${LIB} ${Boost_LIBRARIES})
set(TARGET "KMC")
add_library(${TARGET} INTERFACE)
target_link_libraries(${TARGET} ${LIB})
target_include_directories(${TARGET} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/${TARGET}>
)
