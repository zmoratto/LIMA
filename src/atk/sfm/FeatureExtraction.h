#pragma once
#include <sstream>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

using namespace cv;

struct configParams
{
	int detect_extract_method;
	int showResultsFlag;
	int saveResultsFlag;
};

class FeatureExtraction
{
public:
	//Constructor
	FeatureExtraction();

	//Destructor
	~FeatureExtraction();

	//Member Variables
	configParams featureParams;
	FeatureDetector * detector;
	vector<KeyPoint> key_points;
	DescriptorExtractor * extractor;
	Mat point_descriptors;
	int featureMethod;

	//Methods
	void printWarning(string configFilename);
	void printUsage();
	void readImageFiles(string filename, vector<string>& imageFiles);
	void readConfigurationFile(string filename);
	void setDetectExtract(int type);
	void detectKeyPoints(IplImage* image);
	void extractKeyPoints(IplImage* image);
	void showKeyPoints(IplImage* image);
	void process(IplImage* image);
	void clear();
};
