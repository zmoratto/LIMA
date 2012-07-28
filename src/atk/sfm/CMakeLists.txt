set(SFM_SRCS 	SFM.cpp 
		PoseEstimation.cpp 
		FeatureExtraction.cpp 
		../../../common/Tiling.cpp)

add_library(sfmcommon ${SFM_SRCS})

add_executable(sfm_test sfm_test.cpp)
target_link_libraries(sfm_test sfmcommon
                          ${PCL_LIBRARIES}
                          ${OpenCV_LIBS})

#add_executable(save_kinect save_kinect.cc)
#target_link_libraries(save_kinect ${PCL_IO_LIBRARIES})
#                                  ${PCL_VISUALIZATION_LIBRARIES}
#                                  ${PCL_COMMON_LIBRARIES}
#                                  ${OpenCV_LIBS})

install(TARGETS sfm_test DESTINATION bin)

if (BUILD_TESTS)
  add_subdirectory(tests)
endif()