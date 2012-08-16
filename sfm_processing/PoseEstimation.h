#pragma once
#include <vector>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

//Depth Info For Matching
#define NO_DEPTH 0
#define STEREO_DEPTH 1
#define KINECT_DEPTH 2
#define POINT_CLOUD 3

//Weighting For Matches
#define NO_WEIGHT 0
#define SINGLE_WEIGHT 1
#define DOUBLE_WEIGHT 2

//Point Projection struct for SBA
struct pointProjection{
	cv::Point3f point;
	std::vector<cv::Point2f> pixels;
	std::vector<int> frames;
	std::vector<cv::Point3f> allPoints;
};

class PoseEstimation
{
public:
	PoseEstimation();

	~PoseEstimation();

	//Tile points, matches and descriptors
	std::vector<cv::Point3f> currMatchedPoints, prevMatchedPoints;
	std::vector<cv::Point2f> prevMatchedPixels, currMatchedPixels;
	cv::Mat prevObjectData, currObjectData;
	std::vector<cv::KeyPoint> prevKeypoints;

	//Global Points, matches and descriptors
	std::vector<cv::Point3f> globalCurrMatchedPts, globalPrevMatchedPts;
	std::vector<cv::Point2f> globalPrevMatchedPix, globalCurrMatchedPix;
	std::vector<cv::KeyPoint> globalCurrKeyPoints, globalPrevKeyPoints;
	cv::Mat globalCurrDescriptors, globalPrevDescriptors;

	//Vector with index locations of globalMatchedPix which says which matches are valid
	std::vector<int> validPix;

	//Vector with floats telling which matches have which weights
	std::vector<float> matchWeights;

	//Mask which says which matches are inliers and which are outliers
	cv::Mat outlierMask;


	//Rotation and Translation
	CvMat *prevGlobal_R;
	CvMat *prevGlobal_T;
	CvMat *relative_R;
	CvMat *relative_T;
	CvMat *currGlobal_R;
	CvMat *currGlobal_T;

	//3D Points
	std::vector<float> curr_x;
	std::vector<float> curr_y;
	std::vector<float> curr_z;
	std::vector<float> prev_x;
	std::vector<float> prev_y;
	std::vector<float> prev_z;

	//Configuration Parameters
	float nndrRatio;
	int numNN;
	float zDist;
	float xDist;
	float yDist;
	float nbMatches;
	float detThresh;
	int minMatches;
	int homographyMethod;
	int flannCheck;
	float ransacPix;
	float ransacAccuracy;
	float weightedPortion;
	int matchingMethod;

	//Methods
	int computeRelativeR(int numMatches, int numValidMatches);
	void computeRelativeT(int numMatches, int numWeightedMatches);
	int estimateRelativePose();
	void composePose();
	void find_homography();
	void copyGlobal_R_T();
	void resetPrevPos();
	void allocatePCMemory(int width, int height);
	void nearestNeighborMatching(IplImage* image);
	void collect3DMatches(int width);
	void process(int depthInfo, IplImage* image);
	void clear();
	void removeDuplicates(IplImage* image);
	void removeOutliers(int depthInfo, IplImage* image);


};
