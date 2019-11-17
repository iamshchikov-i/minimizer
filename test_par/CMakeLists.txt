set(target ${MINIMIZER_TESTS_PAR})

file(GLOB hdrs "*.h*")
file(GLOB srcs "*.cpp")

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty")

add_executable(${target} ${srcs} ${hdrs})
target_link_libraries(${target} gtest ${MINIMIZER_LIBRARY})
target_link_libraries(${target} ${MPI_LIBRARIES} )
