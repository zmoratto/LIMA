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
	hessianThresh = 800.0;
	nOctaves = 4;
	nOctaveLayers = 2;
	extended = 0;
	upright = 1;
	octaves = 3;
	octaveLayers = 4;
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

//Open list of image files and save in string vector.
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

//Create new detector/extractor objects.
void FeatureExtraction::setDetectExtract(int type)
{

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

//Display keypoints on current image using small rectangles
void FeatureExtraction::showKeyPoints(IplImage* image)
{
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

//Detect key points and extract descriptors from entire image
void FeatureExtraction::process(IplImage* image)
{
	//Clean Up
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

//Clear FeatureExtraction members
void FeatureExtraction::clear()
{
	key_points.clear();
	point_descriptors.release();
}
