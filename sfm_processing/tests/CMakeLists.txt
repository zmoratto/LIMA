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
add_executable(sfm_test sfm_test.cpp)
target_link_libraries(sfm_test sfm)

#FE_test Executable
add_executable(FE_test FeatureExtraction_test.cpp)
target_link_libraries(FE_test sfm)

#sfmDepth Executable
add_executable(sfmDepth_test depthSfm.cpp ../../stereo_processing/CvStereoBMProcessor.cpp ../../string_util.cc)
target_link_libraries(sfmDepth_test sfm)
target_link_libraries(sfmDepth_test ${Boost_LIBRARIES})
