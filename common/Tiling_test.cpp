#include "Tiling.h"
#include <sstream>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

#include "FeatureExtraction.h"

using namespace std;

int main(int argc, char** argv)
{
	vector<string> inputFiles;
	Tiling im1;
	FeatureExtraction feat;
	Tiling im2;
	string inputFilename;
	string outputFilename;
	string configFilename;
	int frameIndex;

	//Read Command Line Arguments
	if(argc != 4)
	{
		im1.printUsage();
		return (-1);
	}

	configFilename = string(argv[1]);
	inputFilename = string(argv[2]);
	outputFilename = string(argv[3]);

	//Read Tile Configuration File
	im1.readTilingConfigFile(configFilename);

	im2.readTilingConfigFile(configFilename);

	im1.readInputFiles(inputFilename, inputFiles);

	feat.featureParams.detect_extract_method = 0;

	for(frameIndex=0;frameIndex<inputFiles.size()-2;frameIndex++)
	{
		im1.image = cvLoadImage(inputFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);

		im1.createTiles();

		cvSetImageROI(im1.image, im1.tiles[0]);

		feat.process(im1.image, 0);

		cvShowImage("Tile 0 Im1", im1.image);

		cvWaitKey(1000);

		im2.image = cvLoadImage(inputFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);

		im2.createTiles();

		cvSetImageROI(im2.image, im2.tiles[0]);

		feat.process(im2.image, 0);

		cvShowImage("Tile 0 Im2", im2.image);

		cvWaitKey(1000);
	}

	return 0;
}

