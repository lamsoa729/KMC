
include_directories(../extern/catch2)

add_executable(test_kmc test_main.cpp)

target_link_libraries(test_kmc PRIVATE KMC)

add_test(NAME kmc_test COMMAND test_kmc)
