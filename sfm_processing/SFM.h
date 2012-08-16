#ifndef SFM_H
#define SFM_H
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc_c.h>

#include "PoseEstimation.h"
#include "FeatureExtraction.h"
#include "../common/Tiling.h"
#include "../stereo_processing/CvStereoBMProcessor.h"
#include "../string_util.h"

#define NO_DEPTH 0
#define STEREO_DEPTH 1
#define KINECT_DEPTH 2
#define POINT_CLOUD 3

struct SFMParams{
	int firstFrame;
	int lastFrame;
	int frameStep;
	int depthInfo;
	int showResultsFlag;
	int saveResultsFlag;
	string resultsDirName;
	string inputDataName;
	string cameraCalibrationFilename;
	string pointCloudFilename;
	string kinectDepthFilename;
	int poseEstimationType;
	string pointProjFilename;
	string resultImageFilename;
	string stereoFilename;
};

class SFM
{
public:

	//Constructor
	SFM(char* configFile, char* inputFilename);

	//Destructor
	~SFM();

	//Members
	PoseEstimation pose;
	FeatureExtraction feat;
	SFMParams configParams;
	cv::Mat cameraMatrix;
	cv::Mat distCoeffs;
	IplImage *prevImage;
	IplImage *currDepth;
	vector<string> pointFiles;
	Tiling referenceTile;
	int currTile;
	int numTiles;
	int imageWidth;
	int imageHeight;
	vector<string> inputFiles, depthFiles;
	cv::Mat rModelImage;

	//Stereo Members
	CvStereoBMProcessorParameters thisStereoParams;
	cv::Mat rotationMat;
	cv::Mat translationMat;
	CvStereoBMProcessor *thisStereo;
	IplImage* stereoLeft, *stereoRight;
	string stereoLeftList, stereoRightList;
	vector<string> stereoLeftFiles, stereoRightFiles;

	//Point projections for SBA
	std::vector<pointProjection> projectedPoints;

	//Methods
	void setUpSFM(char* inputFilename, IplImage* image);
	void printUsage();
	void printWarning(string& filename);
	void readConfigurationFile(string& configurationFilename);
	void printConfigParams();
	int  readCameraCalibrationFile();
	void restoreDefaultParameters();
	void readImageFilenames(string& inputFile);
	void maskImage(IplImage *image);
	void removeKeyPointsWithoutDepth();
	void processTile(IplImage* image, int frameIndex);
	void preprocessing(IplImage* image, int frameIndex);
	void process(IplImage* image, int frameIndex);
	void update(IplImage* image);
	void showGlobalMatches(IplImage* image1, IplImage* image2, int frameIndex);

	void mapping(string &filename, IplImage* image);
	void depthToPointCloud();
	void readDepthFiles(string &filename);
	void savePointCloud(int iteration, IplImage* image);
	void printCurrentGlobal_R_T();
	int ReadStereoConfigFile(string stereoConfigFilename, CvStereoBMProcessorParameters *thisStereoParams);

	//Point Projection -- For SBA
	int searchPointProj(cv::Point2f find);
	void appendPointProj(cv::Point2f pix, cv::Point3f pt, int loc, int frameIndex);
	void savePointProj(int frameIndex);
	void writePointProj(std::string& pointProjFile);
};

#endif
