#ifndef FEATURE_EXTRACTION_H
#define FEATURE_EXTRACTION_H

#pragma once
#include <sstream>
#include <iostream>
#include <fstream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

class FeatureExtraction
{
public:
	//Constructor
	FeatureExtraction();

	//Destructor
	~FeatureExtraction();

	//Member Variables
	cv::FeatureDetector * detector;
	std::vector<cv::KeyPoint> key_points;
	cv::DescriptorExtractor * extractor;
	cv::Mat point_descriptors;
	int featureMethod;

	//Surf Parameters
	double hessianThresh;
	int nOctaves;
	int nOctaveLayers;
	bool extended;
	bool upright;
	int octaves;
	int octaveLayers;

	//Methods
	void readImageFiles(std::string filename, std::vector<std::string>& imageFiles);
	void setDetectExtract(int type);
	void detectKeyPoints(IplImage* image);
	void extractKeyPoints(IplImage* image);
	void showKeyPoints(IplImage* image);
	void process(IplImage* image);
	void clear();
};

#endif
