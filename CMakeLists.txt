cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ATK)

find_package(PCL 1.5 REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${OPENCV_INCLUDE_DIRS})
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

option(BUILD_TESTS "Build tests" ON)

if (BUILD_TESTS)
  enable_testing()
endif()

add_subdirectory(src/atk/sfm)