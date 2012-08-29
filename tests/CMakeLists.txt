project(lima_tests)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")

set( PROJECT_SRC_DIR ".." )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )

if( CMAKE_COMPILER_IS_GNUCXX )
	set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )

find_package( VisionWorkbench )
if (NOT VW_FOUND)
	message(FATAL_ERROR "Package VisionWorkbench required, but not found.")
endif (NOT VW_FOUND)

find_package( ASP )
if (NOT ASP_FOUND)
	message(FATAL_ERROR "Package ASP required, but not found.")
endif (NOT ASP_FOUND)

find_package( GTest )
if (NOT GTEST_FOUND)
        message(FATAL_ERROR "Package GTEST required, but not found.")
endif (NOT GTEST_FOUND)

## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/..
  ${PROJECT_SOURCE_DIR}
  ${VW_INCLUDE_DIR}
  ${ASP_INCLUDE_DIR}
  ${GTEST_INCLUDE_DIR}
  )
# ISIS include directories
message(STATUS ${GTEST_INCLUDE_DIR})
message(STATUS ${GTEST_LIBRARIES})

set( COMMON_SRCS ../lidar_tracks/tracks.cc ../match.cc ../display.cc ../featuresLOLA.cc ../util.cc )

add_library( lima ${COMMON_SRCS})
target_link_libraries(lima ${VW_LIBRARIES} )
target_link_libraries(lima ${ASP_LIBRARIES} )

# Gets linker error without tracks.cc, even though it should be unneeded...
# I cannot figure out why, but this will do for now
add_executable(test_util Test_util.cc ../util.cc ../lidar_tracks/tracks.cc)
target_link_libraries(test_util ${VW_LIBRARIES} )
target_link_libraries(test_util ${ASP_LIBRARIES} )
target_link_libraries(test_util ${GTEST_LIBRARIES} )

add_executable(test_tracks Test_tracks.cc )
target_link_libraries(test_tracks lima )
target_link_libraries(test_tracks ${VW_LIBRARIES} )
target_link_libraries(test_tracks ${ASP_LIBRARIES} )
target_link_libraries(test_tracks ${GTEST_LIBRARIES} )

add_executable(test_icp Test_icp.cc ../icp.cc )
target_link_libraries(test_icp lima )
target_link_libraries(test_icp ${VW_LIBRARIES} )
target_link_libraries(test_icp ${ASP_LIBRARIES} )
target_link_libraries(test_icp ${GTEST_LIBRARIES} )

