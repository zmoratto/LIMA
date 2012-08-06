#include <sstream>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

#include "../FeatureExtraction.h"
#include <opencv/cv.h>

#include <vector>

using namespace cv;
using namespace std;

void printUsage();
void printWarning(string configFilename);
void readImageFiles(string filename, vector<string>& imageFiles);

int main(int argc, char** argv)
{
	string configFilename;
	string inputFilename;
	string outputDir;
	vector<string> imageFiles;
	FeatureExtraction feat;
	string outputFilename;
	ofstream fout;
	//Command line arguments
	if(argc!=4)
	{
		feat.printUsage();
		return (0);
	}

	feat.featureMethod = atoi(argv[1]);
	inputFilename = string(argv[2]);
	outputDir = string(argv[3]);

	feat.readImageFiles(inputFilename, imageFiles);
	feat.hessianThresh = 800.0;
	feat.nOctaves = 4;
	feat.nOctaveLayers = 2;
	feat.extended = 0;
	feat.upright = 1;
	feat.octaves = 3;
	feat.octaveLayers = 4;

	//Main Loop
	for(int i=0; i<1; i++)
	{
		IplImage * image = cvLoadImage(imageFiles[i].c_str(), CV_LOAD_IMAGE_UNCHANGED);

		//Process Image
		feat.process(image);
		feat.showKeyPoints(image);

		outputFilename = outputDir + imageFiles[i];
		cvSaveImage(outputFilename.c_str(), image);

		cvReleaseImage(&image);
	}

	return 0;
}

