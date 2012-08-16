project(MY_GRAND_PROJECT)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
find_package(PCL 1.5 REQUIRED)
find_package( OpenCV REQUIRED )
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
add_executable(sfmDepth_test depthSfm.cpp ../../stereo_processing/CvStereoBMProcessor.cpp ../../string_util.cc)
target_link_libraries(sfmDepth_test sfm)
target_link_libraries(sfmDepth_test ${Boost_LIBRARIES})


project(mosaic_project)
set( PROJECT_SRC_DIR ".." )
set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules ${CMAKE_MODULE_PATH} )
#possible to edit next line in place of setting environment variable
set(ENV{OALROOT} "~/oal")
find_package( OAL )
if( OAL_FOUND )
   include_directories(${OPENCV_INCLUDE_DIRS} ${OAL_INCLUDE_DIRS})
   add_executable(mosaic ../mosaic.cpp ../../common/opencv_pds.cpp ../../common/pds_read.cpp ../../CAHV-to-pinhole.cc)
   target_link_libraries(mosaic ${OpenCV_LIBS})
   target_link_libraries(mosaic ${OAL})
endif( OAL_FOUND )
if( NOT OAL_FOUND )
   message( "skipping target mosaic" )
endif( NOT OAL_FOUND )


