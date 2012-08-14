project(lidar2img)
cmake_minimum_required(VERSION 2.6)
add_definitions(-arch i386)
add_definitions(-DQT_NO_KEYWORDS)

set( PROJECT_SRC_DIR ".." )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )

if( CMAKE_COMPILER_IS_GNUCXX )
  set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )


find_package(Qt)

find_package( VisionWorkbench )
if (NOT VW_FOUND)
	message(FATAL_ERROR "Package VisionWorkbench required, but not found.")
endif (NOT VW_FOUND)

find_package( ASP )
if (NOT ASP_FOUND)
	message(FATAL_ERROR "Package ASP required, but not found.")
endif (NOT ASP_FOUND)

find_package( Boost REQUIRED system thread filesystem program_options )
if (NOT Boost_FOUND)
	message(FATAL_ERROR "Package Boost required, but not found.")
endif (NOT Boost_FOUND)

## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/
  ${VW_INCLUDE_DIR}
  ${ASP_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  )
# ISIS include directories

set( COMMON_SRCS ../tracks.cc ../match.cc ../coregister.cc ../display.cc ../weights.cc ../featuresLOLA.cc ../util.cc )

add_library( lima ${COMMON_SRCS})
target_link_libraries(lima ${Boost_LIBRARIES} )
target_link_libraries(lima ${VW_LIBRARIES} )
target_link_libraries(lima ${ASP_LIBRARIES} )

add_executable(lidar2img lidar2img.cc )
target_link_libraries(lidar2img lima )
SET_TARGET_PROPERTIES(lima PROPERTIES LINK_FLAGS "-arch i386")
SET_TARGET_PROPERTIES(lidar2img PROPERTIES LINK_FLAGS "-arch i386")
