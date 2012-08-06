######################################################################
# Find script for ASP
#
# Assumes that ASP and ISIS are installed though BinaryBuilder.
# The binary builder base system is assumed to be in the directory
# above ISIS_ROOT.
#
# Input Environment Variables:
# ASPROOT                  : Root of ASP Install
# ISISROOT                 : Root of ISIS install
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
	libQtCore.so.4
	libQtGui.so.4
	libQtNetwork.so.4
	libQtSql.so.4
	libQtSvg.so.4
	libQtXml.so.4
	libQtXmlPatterns.so.4
	libQtWebKit.so.4
	#	boost_filesystem-mt
	#	boost_system-mt
	#	boost_thread-mt
	#	boost_program_options-mt
	superlu
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
	message(ERROR "    ISIS not found. Did you set the ISIS_ROOT environment variable?")
	return()
endif(NOT ISIS_INCLUDE_H)
string(REGEX REPLACE "/[^/]*/[^/]*$" "" ISIS_ROOT_DIR ${ISIS_INCLUDE_H} )

find_file( BASE_INCLUDE_H "include/gdal.h" "$ENV{ISISROOT}/.." NO_DEFAULT_PATH)
if(NOT BASE_INCLUDE_H)
	message(ERROR "    BaseSystem not found. Did you install Isis and set ISIS_ROOT to a directory inside the BinaryBuilder's BaseSystem directory?")
	return()
endif(NOT BASE_INCLUDE_H)
string(REGEX REPLACE "/[^/]*/[^/]*$" "" BASE_SYSTEM_ROOT_DIR ${BASE_INCLUDE_H} )

mark_as_advanced(ASP_INCLUDE_H ISIS_INCLUDE_H BASE_INCLUDE_H)

set( ASP_INCLUDE_DIR 
	${ASP_ROOT_DIR}/include
	${BASE_SYSTEM_ROOT_DIR}/include
	#${BASE_SYSTEM_ROOT_DIR}/include/boost-1_46_1/ 
	${BASE_SYSTEM_ROOT_DIR}/noinstall/include
	${BASE_SYSTEM_ROOT_DIR}/noinstall/include/QtCore
	${ISIS_ROOT_DIR}/3rdParty/include
	${ISIS_ROOT_DIR}/inc
)

set( ASP_LIBRARY_DIR ${ASP_ROOT_DIR}/lib ${ISIS_ROOT_DIR}/lib ${ISIS_ROOT_DIR}/3rdParty/lib ${BASE_SYSTEM_ROOT_DIR}/lib )

LIST(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".so.4") # QT Libraries end with .so.4, need this for find_library to work
foreach(LIB ${ALL_ASP_LIBRARIES})
	set(BLIB BLIB-NOTFOUND) # if we don't do this find_library caches the results
	find_library(BLIB ${LIB} PATHS ${ASP_LIBRARY_DIR} NO_DEFAULT_PATH)
	if(NOT BLIB)
		message(ERROR "    Could not find library ${LIB}.")
		return()
	endif(NOT BLIB)
	set(ASP_LIBRARIES ${ASP_LIBRARIES} ${BLIB} )
endforeach(LIB ${ALL_ASP_LIBRARIES})
LIST(REMOVE_ITEM CMAKE_FIND_LIBRARY_SUFFIXES ".so.4")

set(ASP_FOUND True)

