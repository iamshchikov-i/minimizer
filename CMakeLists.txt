cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME minimizer)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(MINIMIZER_LIBRARY "${PROJECT_NAME}")
set(MINIMIZER_TESTS   "test_${PROJECT_NAME}")
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
add_subdirectory(test)

message( STATUS "")
message( STATUS "General configuration for ${PROJECT_NAME}")
message( STATUS "======================================")
message( STATUS "")
message( STATUS "   Configuration: ${CMAKE_BUILD_TYPE}")
message( STATUS "")

  
