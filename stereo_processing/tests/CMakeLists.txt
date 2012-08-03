project(stereo_project)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
find_package( OpenCV REQUIRED )
find_package( Boost REQUIRED COMPONENTS program_options )
include_directories(
${OPENCV_INCLUDE_DIRS}
${Boost_INCLUDE_DIR}
)

#project(stereo_project)
add_executable(stereo_test CvStereoTest.cpp ../CvStereoProcessor.cpp ../../string_util.cc)
target_link_libraries(stereo_test ${OpenCV_LIBS})
target_link_libraries(stereo_test ${Boost_LIBRARIES})

#project(sgbm_stereo_test)
add_executable(sgbm_stereo_test CvStereoSGBMTest.cpp ../CvStereoSGBMProcessor.cpp ../../common/Tiling.cpp ../../string_util.cc)
target_link_libraries(sgbm_stereo_test ${OpenCV_LIBS})
target_link_libraries(sgbm_stereo_test ${Boost_LIBRARIES})
target_link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY})

#project(bm_stereo_test)
add_executable(bm_stereo_test CvStereoBMTest.cpp ../CvStereoBMProcessor.cpp ../../common/Tiling.cpp ../../string_util.cc)
target_link_libraries(bm_stereo_test ${OpenCV_LIBS})
target_link_libraries(bm_stereo_test ${Boost_LIBRARIES})
