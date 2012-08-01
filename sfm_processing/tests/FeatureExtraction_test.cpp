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

	configFilename = string(argv[1]);
	inputFilename = string(argv[2]);
	outputDir = string(argv[3]);

	feat.readImageFiles(inputFilename, imageFiles);

	//Read configuration file
	feat.readConfigurationFile(configFilename);
	feat.featureMethod = feat.featureParams.detect_extract_method;

	//Main Loop
	for(int i=1; i<2; i++)
	{
		IplImage * image = cvLoadImage(imageFiles[i].c_str(), CV_LOAD_IMAGE_UNCHANGED);

		//Process Image
		feat.process(image);

		if(feat.featureParams.saveResultsFlag)
		{
			outputFilename = outputDir + imageFiles[i];
			cvSaveImage(outputFilename.c_str(), image);
		}

		outputFilename = string("KEYPOINTS_TEST.txt");
		fout.open(outputFilename.c_str());
		for(int j=0;j<feat.key_points.size();j++)
		{
			fout << feat.key_points[j].pt.x << " " << feat.key_points[j].pt.y << endl;
		}
		fout.close();
		cout << "Num KeyPoints: " << feat.key_points.size() << endl;
		cvReleaseImage(&image);
	}

	return 0;
}

