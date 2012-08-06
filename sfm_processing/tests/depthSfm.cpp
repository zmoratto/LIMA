#include <iostream>
#include <time.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

#include "../SFM.h"

using namespace std;

int main(int argc, char** argv)
{
	IplImage *image = NULL;
	IplImage im;
	cv::Mat rModelImage;
	int frameIndex, index;
	ofstream fout;
	string tempString;

	//Read Command Line Arguments
	if(argc != 4)
	{
		cout << endl;
		cout << "******************** USAGE ********************" << endl;
		cout << "./sfm_test <configFile> <inputFile> <outputDir>" << endl;
		cout << "***********************************************" << endl;
		cout << endl;
		return (-1);
	}

	//Set Up SFM
	SFM sfmTest(argv[1], argv[2]);
	if(sfmTest.configParams.depthInfo == STEREO_DEPTH)
	{
		image = cvLoadImage(sfmTest.stereoLeftFiles[0].c_str(), 0);
		sfmTest.setUpSFM(NULL, image);
	}
	else
	{
		image = cvLoadImage(sfmTest.inputFiles[0].c_str(), 0);
		sfmTest.setUpSFM(argv[2], image);
	}

	//Loop through all images in sequence
	for(frameIndex = sfmTest.configParams.firstFrame; frameIndex <= sfmTest.configParams.lastFrame; frameIndex++)
	{
		cout << endl << endl <<"Frame " << frameIndex << endl;

		if(sfmTest.configParams.depthInfo != STEREO_DEPTH)
			image = cvLoadImage(sfmTest.inputFiles[frameIndex].c_str(), 0);

		sfmTest.process(image, frameIndex);
	}

	return 0;
}

