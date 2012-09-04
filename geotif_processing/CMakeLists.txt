project(geotif_processing)
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


## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/..
  ${VW_INCLUDE_DIR}
  )

#set( COMMON_SRCS ../lidar_tracks/tracks.cc ../match.cc ../lidar_tracks/display.cc ../lidar_tracks/featuresLOLA.cc ../util.cc )

add_executable(geotif_resample geotif_resample.cc ../common/string_util.cc)
target_link_libraries(geotif_resample ${VW_LIBRARIES} )



