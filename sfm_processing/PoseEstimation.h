#pragma once
#include <vector>

#include "pcl/common/common_headers.h"
#include <pcl/correspondence.h>
#include <pcl/registration/transformation_estimation_svd.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>

#define NO_DEPTH 0
#define STEREO_DISP 1
#define KINECT_DEPTH 2
#define POINT_CLOUD 3

#define NO_WEIGHT 0
#define SINGLE_WEIGHT 1
#define DOUBLE_WEIGHT 2

using namespace cv;

using namespace std;

struct pointProjection{
	Point3f point;
	vector<Point2f> pixels;
	vector<int> frames;
	vector<Point3f> allPoints;
};

class PoseEstimation
{
public:
	//Constructor
	PoseEstimation();

	//Destructor
	~PoseEstimation();

	//Tile Parameters
	vector<cv::Point3f> currMatchedPoints, prevMatchedPoints;
	vector<cv::Point2f> prevMatchedPixels, currMatchedPixels;
	Mat prevObjectData;
	vector<KeyPoint> prevKeypoints;
	Mat currObjectData;

	//Global Parameters
	vector<cv::Point3f> globalCurrMatchedPts, globalPrevMatchedPts;
	vector<cv::Point2f> globalPrevMatchedPix, globalCurrMatchedPix;
	vector<KeyPoint> globalCurrKeyPoints, globalPrevKeyPoints;
	Mat globalCurrDescriptors, globalPrevDescriptors;
	vector<int> validPix;

	vector<float> matchWeights;

	Mat outlierMask;

	vector<pointProjection> projectedPoints;

	//Rotation and Translation
	CvMat *prevGlobal_R;
	CvMat *prevGlobal_T;
	CvMat *relative_R;
	CvMat *relative_T;
	CvMat *currGlobal_R;
	CvMat *currGlobal_T;

	//3D Points
	vector<float> curr_x;
	vector<float> curr_y;
	vector<float> curr_z;
	vector<float> prev_x;
	vector<float> prev_y;
	vector<float> prev_z;

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

	//Methods
	int estimateRelativePose();
	void composePose();
	void printCurrentGlobal_R_T();
	void clearMatchedPoints();
	void clearMatchedPixels();
	void push_prev3Dpt(int index);
	void push_curr3Dpt(int index);
	int find_curr_index(int loc, int width);
	int find_prev_index(int loc, int width);
	void find_homography();
	void mapping(string &filename, IplImage* image);
	void copyGlobal_R_T();
	void resetPrevPos();
	void depthToPointCloud(IplImage *depthImage, cv::Mat camIntrinsicMatrix);
	void dispToPointCloud(cv::Mat dispMat,  cv::Mat Q);
	void readDepthFiles(string &filename, vector<string>& depthFiles);
	void savePointCloud(vector<string>& filenames, int iteration, IplImage* image);
	void resizeXYZVectors(int width, int height);
	void nearestNeighborMatching();
	void collect3DMatches(int width);
	void process(int depthInfo, IplImage* image);
	void clear();

	//Point Projection -- For SBA
	int searchPointProj(Point2f find);
	void appendPointProj(Point2f pix, Point3f pt, int loc, int frameIndex);
	void savePointProj(int frameIndex, int firstFrame);
	void writePointProj(string& pointProjFile);

	void removeDuplicates(IplImage* image);
};
