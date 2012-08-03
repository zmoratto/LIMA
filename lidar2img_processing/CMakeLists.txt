project(lidar2img)
cmake_minimum_required(VERSION 2.6)

set( PROJECT_SRC_DIR ".." )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )

if( CMAKE_COMPILER_IS_GNUCXX )
  set( WARN_FLAG_ADDITIONS "-fno-strict-aliasing" )
endif( CMAKE_COMPILER_IS_GNUCXX )

find_package( VisionWorkbench )
find_package( ASP )

## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/
  ${VW_INCLUDE_DIR}
  ${ASP_INCLUDE_DIR}
  )
# ISIS include directories

set( COMMON_SRCS ../tracks.cc ../match.cc ../coregister.cc ../display.cc ../weights.cc ../featuresLOLA.cc ../util.cc )

add_library( lima ${COMMON_SRCS})
target_link_libraries(lima ${VW_LIBRARIES} )
target_link_libraries(lima ${ASP_LIBRARIES} )

add_executable(lidar2img ../brute_force_align.cc )
target_link_libraries(lidar2img lima )

