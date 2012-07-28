#include <gtest/gtest.h>

#include <sfm.h>

TEST(TestIO, io_test) {
	//read the config files - START
	SFMParams configParams;
	cv::Mat camIntrinsicMatrix, camDistorsionsCoeffs, Q;

//	EXPECT_EQ(ReadConfigurationFile("", configParams), 0);
//	EXPECT_EQ(ReadCameraCalibrationFile("", camIntrinsicMatrix, camDistorsionsCoeffs), 0);

	/*	
	if ( ReadConfigurationFile(configurationFilename, configParams) ) {
		PrintConfigParams(configParams);
	} else exit(EXIT_FAILURE);
	
	if ( ReadCameraCalibrationFile(cameraCalibrationFilename, camIntrinsicMatrix, camDistorsionsCoeffs) ) {
		cout << "Intrinsic Camera Matrix is " << endl << camIntrinsicMatrix << endl;
		cout << "Camera Distortion Coefficients are " << camDistorsionsCoeffs << "." << endl;
	} else exit(EXIT_FAILURE);
	
	if ( ReadStereoCalibrationFile(stereoCalibrationFilename, Q) ) {
		cout << "Q matrix is " << endl << Q << endl;
	} else exit(EXIT_FAILURE);
	EXPECT_EQ(5, 5);
*/
}
