cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME minimizer)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})

set(MINIMIZER_LIBRARY "${PROJECT_NAME}")
set(MINIMIZER_TESTS_SEQ   "test_${PROJECT_NAME}_seq")
set(MINIMIZER_TESTS_PAR   "test_${PROJECT_NAME}_par")
set(MINIMIZER_TESTS_AGMND  "test_${PROJECT_NAME}_AGMND")
set(MINIMIZER_TESTS_AGP_D  "test_${PROJECT_NAME}_AGP_D")
set(MINIMIZER_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/include")

include_directories("${MINIMIZER_INCLUDE}" gtest)

find_package(MPI)
if(NOT MPI_FOUND)
		message(SEND_ERROR "Failed to find MPI")
		return()
else()
		include_directories(${MPI_INCLUDE_PATH})
endif()

add_subdirectory(src)
add_subdirectory(samples)
add_subdirectory(gtest)
add_subdirectory(test_seq)
add_subdirectory(test_par)
add_subdirectory(test_AGMND)
add_subdirectory(test_AGP_D)


message( STATUS "")
message( STATUS "General configuration for ${PROJECT_NAME}")
message( STATUS "======================================")
message( STATUS "")
message( STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message( STATUS "")

  
