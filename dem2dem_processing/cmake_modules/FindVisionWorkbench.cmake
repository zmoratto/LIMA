######################################################################
# Find script for VisionWorkbench in releases
#
# Input variables:
# ENV variable VWROOT
# ENV variable BASESYSTEMROOT, if desired
#
# Output Variables:
# -----------------
# VW_FOUND                : TRUE if search succeded
# VW_INCLUDE_DIR          : include path
# VW_LIBRARY_DIR          : library path
# VW_LIBRARIES            : All VisionWorkbench libraries that were found
# VW_vwCamera_LIBRARY
# VW_vwCartography_LIBRARY
# VW_vwCore_LIBRARY
# VW_vwFileIO_LIBRARY
# VW_Geometry_LIBRARY
# VW_vwHDR_LIBRARY
# VW_vwImage_LIBRARY
# VW_vwInterestPoint_LIBRARY
# VW_vwMath_LIBRARY
# VW_vwMosaic_LIBRARY
# VW_vwStereo_LIBRARY
# 
######################################################################

set(VW_ROOT $ENV{VWROOT})

include("GetLibraryList")

# Set parent directory as a search location
string(REGEX REPLACE "/[^/]*$" "" PROJ_SRC_PARENT ${PROJECT_SOURCE_DIR})

# these are the modules we'll look for
set( VW_MODULE_NAMES 
  Camera        
  Cartography   
  Core       
  FileIO        
  Geometry          
  HDR           
  Image         
  InterestPoint 
  Math
  Mosaic        
  Stereo        
)

set(     PACKAGE VisionWorkbench )
set( PACKAGE_DIR VisionWorkbench )

find_file( VW_INCLUDE_H "include/vw/vw.h"
  ${VW_ROOT}
)
if(NOT VW_INCLUDE_H)
	message(ERROR "   VW not found. Did you set the VWROOT environment variable?")
	return()
endif(NOT VW_INCLUDE_H)

if (DEFINED ENV{BASESYSTEMROOT})
	SET(BASE_SYSTEM_INCLUDE $ENV{BASESYSTEMROOT}/include)
endif (DEFINED ENV{BASESYSTEMROOT})

# Set the root to 3 directories above vw.h
#-----------------------------------------
string(REGEX REPLACE "/[^/]*/[^/]*/[^/]*$" "" VW_ROOT_DIR ${VW_INCLUDE_H} )


if( VW_ROOT_DIR ) 
  mark_as_advanced(VW_INCLUDE_H)

  set( VW_LIBRARY_NAMES "")
  foreach( MODULE ${VW_MODULE_NAMES} )
    set( VW_LIBRARY_NAMES ${VW_LIBRARY_NAMES} vw${MODULE} )
  endforeach( MODULE ${VW_MODULE_NAMES} )
  
  set( VW_INCLUDE_DIR ${VW_ROOT_DIR}/include ${Boost_INCLUDE_DIR} ${BASE_SYSTEM_INCLUDE})
  set( VW_LIBRARY_DIR ${VW_ROOT_DIR}/lib )
  
  get_library_list(VW ${VW_LIBRARY_DIR} "d" "${VW_LIBRARY_NAMES}")
  get_library_imports(vw "${VW_LIBRARY_DIR}" "${VW_LIBRARY_NAMES}")

  if(WIN32) 
    # create a variable to hint static linkage, following Boost's example
    # this is not currently honored, but we use the varaiable to check
    # whether to add VW_ALL_DYN_LINK to definitions (required for proper
    # dllexport/dllimport on win32)
    set( VW_USE_STATIC_LIBS FALSE CACHE BOOL "hint to link statically (not currently honored)")
    mark_as_advanced( VW_USE_STATIC_LIBS )

    if( NOT VW_USE_STATIC_LIBS )
      add_definitions( -DVW_ALL_DYN_LINK )
      message(STATUS "  Adding VW_ALL_DYN_LINK to compile defines")
    endif( NOT VW_USE_STATIC_LIBS )

    # On windows, although the library dependency is there
    # for transitive linking, the full path is not :(
    link_directories(${Boost_LIBRARY_DIRS})
  endif(WIN32)

  # If we build VW with cmake and we find the import information,
  # we do not need the following 
  if( NOT VW_HAS_IMPORTS )
    # Transitive linking doesn't seem to work very well, on mac
    # so we manually add some extra dependencies
    #if( APPLE ) 
    #  find_library( VW_VECLIB_LIBRARIES Accelerate )
    #  if( VW_VECLIB_LIBRARIES )
    #    set(VW_LIBRARIES ${VW_LIBRARIES} ${VW_VECLIB_LIBRARIES})
    #    set(VW_vwMath_LIBRARY ${VW_vwMath_LIBRARY} ${VW_VECLIB_LIBRARIES})
    #    message( STATUS "  Adding Accelerate framework to vwMath" )
    #    # add boost to vwCore and LIBRARIES
    #    if(NOT Boost_FOUND)
    #      find_package(Boost 1.35.0 COMPONENTS program_options filesystem system thread)
    #    endif(NOT Boost_FOUND)
    #    set(_BOOST_PACKAGES THREAD FILESYSTEM PROGRAM_OPTIONS SYSTEM )
    #    foreach( _BOOST_PACKAGE ${_BOOST_PACKAGES} )
    #      message(STATUS "  Adding Boost_${_BOOST_PACKAGE}_LIBRARY to VW_LIBRARIES and vwCore")
    #      if(Boost_${_BOOST_PACKAGE}_LIBRARY)
    #        set(VW_LIBRARIES      ${VW_LIBRARIES}      ${Boost_${_BOOST_PACKAGE}_LIBRARY} )
    #        set(VW_vwCore_LIBRARY ${VW_vwCore_LIBRARY} ${Boost_${_BOOST_PACKAGE}_LIBRARY} )
    #      else(Boost_${_BOOST_PACKAGE}_LIBRARY)
    #        message(STATUS "   **WARNING** Boost_${_BOOST_PACKAGE}_LIBRARY value is \"${Boost_${_BOOST_PACKAGE}_LIBRARY}\"")
    #        message(STATUS "               but it is needed for VisionWorkbench")
    #      endif(Boost_${_BOOST_PACKAGE}_LIBRARY)
    #    endforeach( _BOOST_PACKAGE ${_BOOST_PACKAGES} )

    #  endif( VW_VECLIB_LIBRARIES )
    #endif( APPLE )
  endif( NOT VW_HAS_IMPORTS )
  
  ## all done
  set( VW_FOUND TRUE )
  
else( VW_ROOT_DIR )
  message(STATUS "  VisionWorkbench NOT found!!! Could not find vw.h")
endif( VW_ROOT_DIR )

if (NOT Boost_FOUND)
	set(ENV{BOOST_ROOT} $ENV{BASESYSTEMROOT}) # use base system root directory
	find_package( Boost REQUIRED filesystem system thread program_options )
	if (NOT Boost_FOUND)
		message(ERROR "Package Boost required for Vision Workbench, but not found.")
		return()
	endif (NOT Boost_FOUND)
	set(VW_INCLUDE_DIR ${VW_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
	set(VW_LIBRARIES ${VW_LIBRARIES} ${Boost_LIBRARIES})
endif (NOT Boost_FOUND)

message(STATUS "  ${Boost_INCLUDE_DIR}")
message(STATUS "  ${Boost_LIBRARIES}")