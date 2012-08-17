######################################################################
# Find script for ASP
#
# Assumes that ASP and ISIS are installed though BinaryBuilder.
# The binary builder base system is assumed to be in the directory
# above ISIS_ROOT. This package also includes the proper QT version
# and boost version for ASP, isis3, and superlu if it is available.
#
# Input Environment Variables:
# ASPROOT                  : Root of ASP Install
# ISISROOT                 : Root of ISIS install
# BASESYSTEMROOT           : Root of BaseSystem, if ISISROOT isn't a subdirectory
#
# Output Variables:
# -----------------
# ASP_FOUND                : TRUE if search succeded
# ASP_INCLUDE_DIR          : include path
# ASP_LIBRARIES            : All ASP libraries that were found
# ASP_ROOT_DIR             : Root directory of ASP Install
# BASE_SYSTEM_ROOT_DIR     : Root directory of Base System Install
# ISIS_ROOT_DIR            : Root directory of Isis Install
# 
######################################################################
set(     PACKAGE ASP )
set( PACKAGE_DIR ASP )

set(ALL_ASP_LIBRARIES
	aspIsisIO
	isis3
)

set(OPTIONAL_ASP_LIBRARIES
	superlu # on mac this is not needed?
)

set(ARBITRARY_QT_LIBRARIES
	libQtCore.so.4
	libQtGui.so.4
	libQtNetwork.so.4
	libQtSql.so.4
	libQtSvg.so.4
	libQtXml.so.4
	libQtXmlPatterns.so.4
	libQtWebKit.so.4
)

set(ARBITRARY_BOOST_LIBRARIES
	libboost_system-mt-1_46_1.so
	libboost_filesystem-mt-1_46_1.so
	libboost_thread-mt-1_46_1.so
	libboost_program_options-mt-1_46_1.so
)

# Set the root to 3 directories above Core.h
# Look for files to confirm that paths are correct
find_file( ASP_INCLUDE_H "include/asp/Core.h" $ENV{ASPROOT} NO_DEFAULT_PATH)
if(NOT ASP_INCLUDE_H)
	message(ERROR "   ASP not found. Did you set the ASPROOT environment variable?")
	return()
endif(NOT ASP_INCLUDE_H)
string(REGEX REPLACE "/[^/]*/[^/]*/[^/]*$" "" ASP_ROOT_DIR ${ASP_INCLUDE_H} )

find_file( ISIS_INCLUDE_H "inc/Isis.h" $ENV{ISISROOT} NO_DEFAULT_PATH)
if(NOT ISIS_INCLUDE_H)
	message(ERROR "    ISIS not found. Did you set the ISISROOT environment variable?")
	return()
endif(NOT ISIS_INCLUDE_H)
string(REGEX REPLACE "/[^/]*/[^/]*$" "" ISIS_ROOT_DIR ${ISIS_INCLUDE_H} )

if (DEFINED ENV{BASESYSTEMROOT})
	set(BASESYSTEMROOT "$ENV{BASESYSTEMROOT}")
else (DEFINED ENV{BASESYSTEMROOT})
	set(BASESYSTEMROOT "$ENV{ISISROOT}/..")
endif (DEFINED ENV{BASESYSTEMROOT})

find_file( BASE_INCLUDE_H "include/gdal.h" ${BASESYSTEMROOT} NO_DEFAULT_PATH)
if(NOT BASE_INCLUDE_H)
	message(ERROR "    BaseSystem not found. Did you install Isis and set ISIS_ROOT to a directory inside the BinaryBuilder's BaseSystem directory?")
	return()
endif(NOT BASE_INCLUDE_H)
string(REGEX REPLACE "/[^/]*/[^/]*$" "" BASE_SYSTEM_ROOT_DIR ${BASE_INCLUDE_H} )

mark_as_advanced(ASP_INCLUDE_H ISIS_INCLUDE_H BASE_INCLUDE_H)

set( ASP_INCLUDE_DIR 
	${ASP_ROOT_DIR}/include
	${BASE_SYSTEM_ROOT_DIR}/include
	${BASE_SYSTEM_ROOT_DIR}/noinstall/include
	${BASE_SYSTEM_ROOT_DIR}/noinstall/include/QtCore
	${ISIS_ROOT_DIR}/3rdParty/include
	${ISIS_ROOT_DIR}/inc
	/opt/local/include # for OS X
)

set( ASP_LIBRARY_DIR ${ASP_ROOT_DIR}/lib ${ISIS_ROOT_DIR}/lib ${ISIS_ROOT_DIR}/3rdParty/lib ${BASE_SYSTEM_ROOT_DIR}/lib )

foreach(LIB ${ALL_ASP_LIBRARIES})
	set(BLIB BLIB-NOTFOUND) # if we don't do this find_library caches the results
	find_library(BLIB ${LIB} PATHS ${ASP_LIBRARY_DIR} NO_DEFAULT_PATH)
	if(NOT BLIB)
		message(ERROR "    Could not find library ${LIB}.")
		return()
	endif(NOT BLIB)
	set(ASP_LIBRARIES ${ASP_LIBRARIES} ${BLIB} )
endforeach(LIB ${ALL_ASP_LIBRARIES})

foreach(LIB ${OPTIONAL_ASP_LIBRARIES})
	set(BLIB BLIB-NOTFOUND) # if we don't do this find_library caches the results
	find_library(BLIB ${LIB} PATHS ${ASP_LIBRARY_DIR} NO_DEFAULT_PATH)
	if(NOT BLIB)
		message(WARNING "    Could not find library ${LIB}.")
	else(NOT BLIB)
		set(ASP_LIBRARIES ${ASP_LIBRARIES} ${BLIB} )
	endif(NOT BLIB)
endforeach(LIB ${OPTIONAL_ASP_LIBRARIES})

LIST(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".so.4") # QT Libraries end with .so.4, need this for find_library to work
set( ARBITRARY_QT_FOUND True )
foreach(LIB ${ARBITRARY_QT_LIBRARIES})
	set(BLIB BLIB-NOTFOUND) # if we don't do this find_library caches the results
	find_library(BLIB ${LIB} PATHS ${ASP_LIBRARY_DIR} NO_DEFAULT_PATH)
	if(NOT BLIB)
		message(WARNING "    Could not find QT library ${LIB}.")
		set( ARBITRARY_QT_FOUND False )
	else(NOT BLIB)
		set(ASP_LIBRARIES ${ASP_LIBRARIES} ${BLIB} )
	endif(NOT BLIB)
endforeach(LIB ${ARBITRARY_QT_LIBRARIES})
LIST(REMOVE_ITEM CMAKE_FIND_LIBRARY_SUFFIXES ".so.4")

if (NOT ARBITRARY_QT_FOUND)
	message(STATUS "    Using system Qt version.")
	find_package(Qt)
	if (NOT Qt_FOUND)
		message(ERROR "Package Qt required for ASP, but not found.")
		return()
	endif (NOT Qt_FOUND)
	set(ASP_INCLUDE_DIR ${ASP_INCLUDE_DIR} ${Qt_INCLUDE_DIR})
	set(ASP_LIBRARIES ${ASP_LIBRARIES} ${Qt_LIBRARIES})
endif (NOT ARBITRARY_QT_FOUND)

set( ARBITRARY_BOOST_FOUND True )
foreach(LIB ${ARBITRARY_BOOST_LIBRARIES})
	set(BLIB BLIB-NOTFOUND) # if we don't do this find_library caches the results
	find_library(BLIB ${LIB} PATHS ${ASP_LIBRARY_DIR} NO_DEFAULT_PATH)
	if(NOT BLIB)
		message(WARNING "    Could not find Boost library ${LIB}.")
		set( ARBITRARY_BOOST_FOUND False )
	else(NOT BLIB)
		set(ASP_LIBRARIES ${ASP_LIBRARIES} ${BLIB} )
	endif(NOT BLIB)
endforeach(LIB ${ARBITRARY_BOOST_LIBRARIES})

if (NOT ARBITRARY_BOOST_FOUND)
	find_package( Boost REQUIRED filesystem system thread program_options )
	if (NOT Boost_FOUND)
		message(ERROR "Package Boost required for ASP, but not found.")
		return()
	endif (NOT Boost_FOUND)
	set(ASP_INCLUDE_DIR ${ASP_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
	set(ASP_LIBRARIES ${ASP_LIBRARIES} ${Boost_LIBRARIES})
endif (NOT ARBITRARY_BOOST_FOUND)

add_definitions(-DQT_NO_KEYWORDS) #avoids conflict if using Boost 1.48

set(ASP_FOUND True)

