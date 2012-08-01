#include "Tiling.h"
#include <sstream>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;

int main(int argc, char** argv)
{
	vector<string> inputFiles;
	Tiling im1;
	Tiling im2;
	string inputFilename;
	string outputFilename;
	string configFilename;
	int frameIndex;
	IplImage* image;

	//Read Command Line Arguments
	if(argc != 3)
	{
		im1.printUsage();
		return (-1);
	}

	configFilename = string(argv[1]);
	inputFilename = string(argv[2]);

	//Read Tile Configuration File
	im1.readTilingConfigFile(configFilename);

	im1.readInputFiles(inputFilename, inputFiles);

	for(frameIndex=0;frameIndex<inputFiles.size()-1;frameIndex++)
	{
		image = cvLoadImage(inputFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);

		im1.setUpReferenceTiles(image);

		im1.process(image);

		for(int i=0; i<im1.tiles.size(); i++)
		{
			cout << "BB " << i << ": [" << im1.tiles[i].x << ", " << im1.tiles[i].y <<  ", " << im1.tiles[i].width << ", " << im1.tiles[i].height << "]" << endl; 
			cvShowImage("Tile", im1.tileImages[i]);
			cvWaitKey(1000);
		}
	}

	return 0;
}

