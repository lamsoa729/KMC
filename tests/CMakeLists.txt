find_package(Boost REQUIRED)
include_directories(catch2 ${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})
add_executable(test_kmc kmc_test.cpp
                        ${CMAKE_SOURCE_DIR}/src/kmc/kmc.cpp
                        ${CMAKE_SOURCE_DIR}/src/kmc/lookup_table.cpp
                        ${CMAKE_SOURCE_DIR}/src/kmc/integrals.cpp
)
target_link_libraries(test_kmc ${LIBS})
add_executable(test_lookup lookup_test.cpp
                        ${CMAKE_SOURCE_DIR}/src/kmc/lookup_table.cpp
                        ${CMAKE_SOURCE_DIR}/src/kmc/integrals.cpp
)
target_link_libraries(test_lookup ${LIBS})
add_test(kmc test_kmc)
add_test(lookup_table test_lookup)
