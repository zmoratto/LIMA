######################################################################
# Find script for SBA
#
# Assumes that environment variable SBAROOT is set to the root of the oal install 
#
# Input Environment Variables:
# SBAROOT                  : Root of SBA Install
#
# Output Variables:
# -----------------
# SBA_FOUND                : TRUE if search succeded
# SBA_INCLUDE_DIRS         : include paths for oal
# SBA                      : library if found
# 
######################################################################

set(SBA_ROOT_DIR $ENV{SBAROOT})
set(SBA_SRC_DIR "${SBA_ROOT_DIR}")
set(SBA_DIR "${SBA_SRC_DIR}")
set(SBA_LIB_PATH  "${SBA_SRC_DIR}")

set(SBA_INCLUDE_DIRS ${SBA_DIR} ${SBA_SRC_DIR} ${SBA_ROOT_DIR})

set(SBA SBA-NOTFOUND) # if we don't do this find_library caches the results
find_library(SBA ${SBA_LIB} PATHS ${SBA_LIB_PATH} ${SBA_ROOT_DIR} )

#if not found search in the sba root directory for libsba.a 
if(NOT SBA)
   find_library(SBA sba PATHS ${SBA_ROOT_DIR})
endif(NOT SBA)

#if still not found print error message
if(NOT SBA)
   message("WARNING: Could not find library ${SBA_LIB}\nDid you set SBAROOT environment variable?")
   set(SBA_FOUND false)
endif(NOT SBA)
if(SBA)
   set(SBA_FOUND true)
endif(SBA)



