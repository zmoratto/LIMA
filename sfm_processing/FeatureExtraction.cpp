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
		delete detector;
	if(extractor != NULL)
		delete extractor;
}

void FeatureExtraction::printWarning(string configFilename)
{
	cout << endl;
	cout << "**************************************************************" << endl;
	cout << "WARNING: Unable to open file " << configFilename << "..." << endl;
	cout << "Using default parameters..." << endl;
	cout << "Feature Detection and Extraction Method: Surf" << endl;
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
		featureParams.detect_extract_method = 0;
		featureParams.showResultsFlag = 1;
		featureParams.saveResultsFlag = 1;
	}
}

void FeatureExtraction::setDetectExtract(int type)
{
	switch(type)
	{
		case 1: //Orb
			detector = new OrbFeatureDetector();
			extractor = new OrbDescriptorExtractor();
			break;
		case 2: //Sift
			detector = new SiftFeatureDetector();
			extractor = new SiftDescriptorExtractor();
			break;
		default: //Surf
			detector = new SurfFeatureDetector(hessianThresh, octaves, octaveLayers, upright);
			extractor = new SurfDescriptorExtractor(nOctaves, nOctaveLayers, extended, upright);
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
	CvPoint pt1, pt2; 
	for (int i = 0; i < key_points.size(); i++)
	{
		pt1.x = key_points[i].pt.x-1;
		pt2.x = key_points[i].pt.x+1; 
		pt1.y = key_points[i].pt.y-1;
		pt2.y = key_points[i].pt.y+1;
		cvRectangle(image, pt1, pt2, CV_RGB(100, 0, 255), 3, 8, 0 );
	}
}

void FeatureExtraction::process(IplImage* image)
{
	vector<KeyPoint> temp;
	int cnt = 0;
	int xLoc, yLoc, width = image->width;
	temp.clear();
	key_points.clear();
	point_descriptors.release();
	if(detector != NULL)
		delete detector;
	if(extractor != NULL)
		delete extractor;
	setDetectExtract(featureMethod);
	detectKeyPoints(image);

	extractKeyPoints(image);
}

void FeatureExtraction::clear()
{
	key_points.clear();
	point_descriptors.release();
}
