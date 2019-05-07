find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})
add_executable(kmc  main.cpp
                    kmc.cpp
                    kmc.hpp
                    integrals.cpp
                    integrals.hpp
                    lookup_table.cpp
                    lookup_table.hpp
)
target_link_libraries(kmc ${LIBS})

