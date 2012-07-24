cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(STEREO_PROCESSING)
find_package( OpenCV REQUIRED )
find_package( Boost REQUIRED)
include_directories(
${OPENCV_INCLUDE_DIRS}
${Boost_INCLUDE_DIR}
)
project(bm_stereo_project)
add_executable(bm_stereo_test CvStereoTest.cpp CvStereoProcessor.cpp)
target_link_libraries(bm_stereo_test ${OpenCV_LIBS})
target_link_libraries(bm_stereo_test ${Boost_LIBS})
project(bm_stereo_project)
add_executable(sbgm_stereo_test CvSBGMStereoTest.cpp CvSBGMStereoProcessor.cpp)
target_link_libraries(sbgm_stereo_test ${OpenCV_LIBS})
target_link_libraries(sbgm_stereo_test ${Boost_LIBS})
