######################################################################
# Find script for OAL
#
# Assumes that environment variable OALROOT is set to the root of the oal install 
#
# Input Environment Variables:
# OALROOT                  : Root of ASP Install
#
# Output Variables:
# -----------------
# OAL_FOUND                : TRUE if search succeded
# ASP_INCLUDE_DIR          : include path
# ASP_LIBRARIES            : All ASP libraries that were found
# ASP_ROOT_DIR             : Root directory of ASP Install
# BASE_SYSTEM_ROOT_DIR     : Root directory of Base System Install
# ISIS_ROOT_DIR            : Root directory of Isis Install
# 
######################################################################
set( PACKAGE OAL )

#set(ENV{OALROOT} "~/oal")

set(OAL_ROOT_DIR $ENV{OALROOT})
set(OAL_SRC_DIR "${OAL_ROOT_DIR}/src")
set(OAL_DIR "${OAL_SRC_DIR}/oal")
set(LABLIB3_DIR "${OAL_SRC_DIR}/lablib3")
set(OAL_LIB_PATH  "${OAL_SRC_DIR}/lib")
set(OAL_LIB "oal_others_c")

set(OAL_INCLUDE_DIRS ${OAL_DIR} ${LABLIB3_DIR} ${OAL_SRC_DIR} ${OAL_ROOT_DIR})

set(OAL OAL-NOTFOUND) # if we don't do this find_library caches the results
find_library(OAL ${OAL_LIB} PATHS ${OAL_LIB_PATH} )

#if not found search in the oal root directory for liboal.a (used in old versions of oal)
if(NOT OAL)
   find_library(OAL oal PATHS ${OAL_ROOT_DIR})
endif(NOT OAL)

#if still not found print error message
if(NOT OAL)
   message(WARNING "Could not find library ${OAL_LIB}\nDid you set OALROOT environment variable")
   message("Cannot make pds_test_read without oal")
   set(OAL_FOUND false)
endif(NOT OAL)
if(OAL)
   set(OAL_FOUND true)
endif(OAL)



