project(MY_GRAND_PROJECT)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
find_package(PCL 1.5 REQUIRED)
find_package(OpenCV REQUIRED )



include_directories(${OPENCV_INCLUDE_DIRS})
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
add_definitions(-g)

#sfm Library
add_library(sfm ../SFM.cpp ../PoseEstimation.cpp ../FeatureExtraction.cpp ../../common/Tiling.cpp)
target_link_libraries(sfm ${PCL_LIBRARIES})
target_link_libraries(sfm ${OpenCV_LIBS})

#sfm_test Executable
#add_executable(sfm_test sfm_test.cpp)
#target_link_libraries(sfm_test sfm)

#FE_test Executable
add_executable(FE_test FeatureExtraction_test.cpp)
target_link_libraries(FE_test sfm)

#sfmDepth Executable
add_executable(sfmDepth_test depthSfm.cpp ../../stereo_processing/CvStereoBMProcessor.cpp ../../common/string_util.cc)
target_link_libraries(sfmDepth_test sfm)
target_link_libraries(sfmDepth_test ${Boost_LIBRARIES})

#eucsbadeo executable
set( PROJECT_SRC_DIR ".." )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )
find_package(SBA)
include(CMakeFindFrameworks)
CMAKE_FIND_FRAMEWORKS(vecLib)
#INCLUDE_DIRECTORIES(..) 
INCLUDE_DIRECTORIES(${SBA_INCLUDE_DIRS}) 
LINK_DIRECTORIES(${SBA_LIB_PATH} ${LAPACKBLAS_DIR})
ADD_EXECUTABLE(eucsbademo eucsbademo.c imgproj.c readparams.c eucsbademo.h readparams.h)
# libraries the demo depends on
IF(HAVE_F2C)  
 TARGET_LINK_LIBRARIES(eucsbademo sba ${LAPACK_LIB} ${BLAS_LIB} ${F2C_LIB})
ELSE(HAVE_F2C)
  TARGET_LINK_LIBRARIES(eucsbademo sba ${LAPACK_LIB} ${BLAS_LIB} ${vecLib_FRAMEWORKS} )
ENDIF(HAVE_F2C)
# make sure that the library is built before the demo
ADD_DEPENDENCIES(eucsbademo sba)

project(mosaic_project)
set( PROJECT_SRC_DIR ".." )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )
#possible to edit next line in place of setting environment variable
#set(ENV{OALROOT} "~/oal")
find_package( OAL )
if( OAL_FOUND )
   include_directories(${OPENCV_INCLUDE_DIRS} ${OAL_INCLUDE_DIRS})
   add_executable(mosaic_test mosaicTest.cpp ../mosaic.cpp ../../common/opencv_pds.cpp ../../common/pds_read.cpp ../../common/string_util.cc ../../camera_models/CAHV-to-pinhole.cc)
   target_link_libraries(mosaic_test ${OpenCV_LIBS})
   target_link_libraries(mosaic_test ${OAL})
endif( OAL_FOUND )
if( NOT OAL_FOUND )
   message( "Skipping target mosaic_test" )
endif( NOT OAL_FOUND )


