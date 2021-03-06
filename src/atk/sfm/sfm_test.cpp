#include <iostream>
#include <time.h>
#include "PoseEstimation.h"
#include "FeatureExtraction.h"
#include "SFM.h"

using namespace std;

int main (int argc, char** argv)
{
	IplImage* image = NULL;
	clock_t t1, t2;
	KeyPoint tilePt;
	int oldRow;

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

	//Read Config Files and set up tiles
	SFM sfmTest(argv[1], argv[2]);

	//Process Each Image Frame
	for(int frameIndex=sfmTest.configParams.firstFrame; frameIndex<=sfmTest.configParams.lastFrame; frameIndex++)
	{
		cout << "******************************************" << endl << endl;
		cout << "FRAME NUM: " << frameIndex << endl;

		image = cvLoadImage(sfmTest.inputFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);
		t1 = clock();
		sfmTest.process(image, frameIndex);
		t2 = clock();

		cout << "Time: " << (double(t2)-double(t1))/CLOCKS_PER_SEC << endl;
		cout << endl << "******************************************" << endl << endl;
	}

	return (0);
}

