set(SFM_SRCS PairwisePose.cpp 
             PoseComposition.cpp 
             display.cpp 
             KinectCalibrate.cpp 
             mapping.cpp 
             io.cpp)

add_library(sfmcommon ${SFM_SRCS})

add_executable(sfm sfm.cpp)
target_link_libraries(sfm sfmcommon
                          ${PCL_LIBRARIES}
                          ${OpenCV_LIBS})

add_executable(pc_vis pc_vis.cpp)
target_link_libraries(pc_vis ${PCL_LIBRARIES})

#add_executable(save_kinect save_kinect.cc)
#target_link_libraries(save_kinect ${PCL_IO_LIBRARIES})
#                                  ${PCL_VISUALIZATION_LIBRARIES}
#                                  ${PCL_COMMON_LIBRARIES}
#                                  ${OpenCV_LIBS})

add_executable(calibration DoKinectCalibration.cpp)
target_link_libraries(calibration ${PCL_LIBRARIES}
                                  ${OpenCV_LIBS})

add_executable(imagelist_creator imagelist_creator.cpp)
target_link_libraries(imagelist_creator ${OpenCV_LIBS})

install(TARGETS sfm pc_vis DESTINATION bin)

if (BUILD_TESTS)
  add_subdirectory(test)
endif()
