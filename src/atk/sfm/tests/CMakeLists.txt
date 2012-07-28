include_directories(${CMAKE_SOURCE_DIR}/thirdparty/gtest)
include_directories(${CMAKE_SOURCE_DIR}/thirdparty/gtest/include)
add_library(gtest ${CMAKE_SOURCE_DIR}/thirdparty/gtest/src/gtest-all.cc)

add_library(gtest_main ${CMAKE_SOURCE_DIR}/thirdparty/gtest/src/gtest_main.cc)
target_link_libraries(gtest_main gtest)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

macro(add_sfm_test test_name)
  add_executable(${test_name} ${test_name}.cxx)
  target_link_libraries(${test_name} sfmcommon
                                     gtest_main
                                     ${PCL_LIBRARIES}
                                     ${OpenCV_LIBS})
  add_test(${test_name} ${test_name})
endmacro()

add_sfm_test(TestIO)
add_sfm_test(TestExample)
add_sfm_test(TestSomethingElse)
