#include <iostream>
#include <time.h>
#include "../SFM.h"

using namespace std;

int main (int argc, char** argv)
{
	IplImage* image = NULL;

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
	SFM sfmTest(argv[1]);
	sfmTest.setUpSFM(argv[2], image);

	//Process Each Image Frame
	for(int frameIndex=sfmTest.configParams.firstFrame; frameIndex<=sfmTest.configParams.lastFrame; frameIndex+=sfmTest.configParams.frameStep)
	{
		cout << "******************************************" << endl << endl;
		cout << "FRAME NUM: " << frameIndex << endl;

		image = cvLoadImage(sfmTest.inputFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);
		sfmTest.process(image, frameIndex);

		cout << endl << "******************************************" << endl << endl;
	}

	//Write point projections for SBA
	sfmTest.pose.writePointProj(sfmTest.configParams.pointProjFilename);

	return (0);
}

