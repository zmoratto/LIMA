project(lidar2img)
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

#find_package( ASP )
#if (NOT ASP_FOUND)
#	message(FATAL_ERROR "Package ASP required, but not found.")
#endif (NOT ASP_FOUND)

## We must add the in-source and out-of-source
## include paths 
##--------------------------------------------
include_directories(
  ${PROJECT_SOURCE_DIR}/..
  ${VW_INCLUDE_DIR}
  #  ${ASP_INCLUDE_DIR}
  )

set( COMMON_SRCS ../icp.cc ../common/string_util.cc pyr_tiling.cc assembler.cc )

add_executable( assembler ${COMMON_SRCS})
target_link_libraries(assembler ${VW_LIBRARIES} )
