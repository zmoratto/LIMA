#include <iostream>
#include "FeatureExtraction.h"
#include <time.h>
#include <vector>

using namespace std;


//Constructor
FeatureExtraction::FeatureExtraction()
{
	detector = NULL;
	extractor = NULL;
	featureMethod = 0;
}

//Destructor
FeatureExtraction::~FeatureExtraction()
{
	if(detector != NULL)
	{
		delete detector;
		detector = NULL;
	}
	if(extractor != NULL)
	{
		delete extractor;
		extractor = NULL;
	}
}

void FeatureExtraction::printWarning(string configFilename)
{
	//Warning for default parameters
	cout << endl;
	cout << "**************************************************************" << endl;
	cout << "WARNING: Unable to open file " << configFilename << "..." << endl;
	cout << "Using default parameters..." << endl;
	cout << "**************************************************************" << endl;
	cout << endl;
}

void FeatureExtraction::printUsage()
{
	cout << endl;
	cout << "********************* USAGE ********************" << endl;
	cout << "./sfm_test <configFile> <inputFile> <outputDir> " << endl;
	cout << "************************************************" << endl;
	cout << endl;
}

void FeatureExtraction::readImageFiles(string filename, vector<string>& imageFiles)
{
	//Open list of image files and save in string vector.
	ifstream fin;
	string temp;
	fin.open(filename.c_str());
	imageFiles.clear();
	while(fin.good())
	{
		getline(fin, temp);
		imageFiles.push_back(temp);
	}

	fin.close();
}

void FeatureExtraction::readConfigurationFile(string filename)
{
	ifstream fin;
	int detect_extract_method;
	int showResultsFlag;
	int saveResultsFlag;
	string identifier;

	//Read configuration file for FE_test
	fin.open(filename.c_str());

	if(fin.is_open())
	{
		fin >> identifier >> detect_extract_method;
		fin >> identifier >> showResultsFlag;
		fin >> identifier >> saveResultsFlag;

		fin.close();

		featureParams.detect_extract_method = detect_extract_method;
		featureParams.showResultsFlag = showResultsFlag;
		featureParams.saveResultsFlag = saveResultsFlag;
	}
	else //Default Parameters
	{
		printWarning(filename);
		featureParams.detect_extract_method = 0; //SURF
		featureParams.showResultsFlag = 1;
		featureParams.saveResultsFlag = 1;
	}
}

void FeatureExtraction::setDetectExtract(int type)
{
	//Create new detector/extractor objects.
	switch(type)
	{
		case 1: //Orb
			detector = new cv::OrbFeatureDetector();
			extractor = new cv::OrbDescriptorExtractor();
			break;
		case 2: //Sift
			detector = new cv::SiftFeatureDetector();
			extractor = new cv::SiftDescriptorExtractor();
			break;
		default: //Surf
			detector = new cv::SurfFeatureDetector(hessianThresh, octaves, octaveLayers, upright);
			extractor = new cv::SurfDescriptorExtractor(nOctaves, nOctaveLayers, extended, upright);
	}
}

//Detect Key Points
void FeatureExtraction::detectKeyPoints(IplImage* image)
{
	detector->detect(image, key_points);
}

//Extract Descriptors
void FeatureExtraction::extractKeyPoints(IplImage* image)
{
	extractor->compute(image, key_points, point_descriptors);
}


//Show Key Points
void FeatureExtraction::showKeyPoints(IplImage* image)
{
	//Display keypoints on current image using small rectangles
	CvPoint pt1, pt2; 
	for (int i = 0; i < key_points.size(); i++)
	{
		pt1.x = key_points[i].pt.x-1;
		pt2.x = key_points[i].pt.x+1; 
		pt1.y = key_points[i].pt.y-1;
		pt2.y = key_points[i].pt.y+1;
		cvRectangle(image, pt1, pt2, CV_RGB(255, 255, 255), 3, 8, 0 );
	}
}

void FeatureExtraction::process(IplImage* image)
{
	//Clean Up
	vector<cv::KeyPoint> temp;
	int cnt = 0;
	int xLoc, yLoc, width = image->width;
	temp.clear();
	key_points.clear();
	point_descriptors.release();
	if(detector != NULL)
	{
		delete detector;
		detector = NULL;
	}
	if(extractor != NULL)
	{
		delete extractor;
		detector = NULL;
	}

	//Detect and extract keypoints from image
	setDetectExtract(featureMethod);
	detectKeyPoints(image);
	extractKeyPoints(image);
}

void FeatureExtraction::clear()
{
	key_points.clear();
	point_descriptors.release();
}
