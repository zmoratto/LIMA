project(lidar2dem)
cmake_minimum_required(VERSION 2.6)

#add_definitions(-DQT_NO_KEYWORDS) #avoids conflict if using Boost 1.48

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

#find_package( Boost REQUIRED system thread filesystem program_options )
#if (NOT Boost_FOUND)
#	message(FATAL_ERROR "Package Boost required, but not found.")
#endif (NOT Boost_FOUND)

## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/..
  ${VW_INCLUDE_DIR}
  ${ASP_INCLUDE_DIR}
  #${Boost_INCLUDE_DIR}
  )
# ISIS include directories

set( COMMON_SRCS ../icp.cc  ../lidar_tracks/tracks.cc  ../util.cc lidar2dem_config.cc lidar2dem.cc )

add_executable( lidar2dem ${COMMON_SRCS})
#target_link_libraries(lidar2dem ${Boost_LIBRARIES} )
target_link_libraries(lidar2dem ${VW_LIBRARIES} )
target_link_libraries(lidar2dem ${ASP_LIBRARIES} )

#IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") # only 32-bit ISIS is available on Mac
#	add_definitions(-arch i386)
#	SET_TARGET_PROPERTIES(lidar2dem PROPERTIES LINK_FLAGS "-arch i386")
#	SET_TARGET_PROPERTIES(lidar2dem PROPERTIES LINK_FLAGS "-arch i386")
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

