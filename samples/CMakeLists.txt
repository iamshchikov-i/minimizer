add_executable(sample main.cpp)
target_link_libraries(sample ${MINIMIZER_LIBRARY})
target_link_libraries(sample ${MPI_LIBRARIES} )
