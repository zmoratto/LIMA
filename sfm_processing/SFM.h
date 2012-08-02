#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/registration/icp.h>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/imgproc/imgproc_c.h>

#include "PoseEstimation.h"
#include "FeatureExtraction.h"
#include "../common/Tiling.h"

#define NO_DEPTH 0
#define KINECT_DEPTH 2
#define POINT_CLOUD 3

using namespace cv;
using namespace std;

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
	int featureMethod;
	int poseEstimationType;
	string pointProjFilename;
	string resultImageFilename;
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
	Mat cameraMatrix;
	Mat distCoeffs;
	Mat Q;
	IplImage *prevImage;
	IplImage *prevDepth, *currDepth;
	string outputDir;
	int sfmInit;
	vector<string> pointFiles;
	Tiling referenceTile;
	Tiling matchingTile;
	int currTile;
	int numTiles;
	int imageWidth;
	int imageHeight;
	vector<string> inputFiles, depthFiles;

	//Methods
	void printUsage();
	void printWarning(string& filename);
	void readConfigurationFile(string& configurationFilename);
	void printConfigParams();
	int  readCameraCalibrationFile();
	void restoreDefaultParameters();
	void readImageFilenames(string& inputFile);
	void processTile(IplImage* image, int frameIndex);
	void preprocessing(IplImage* image, int frameIndex);
	void process(IplImage* image, int frameIndex);
	void update(IplImage* image);
	void clear();
	void showGlobalMatches(IplImage* image1, IplImage* image2, int frameIndex);
	void cleanUp(IplImage*& im1, IplImage*& im2);

	void mapping(string &filename, IplImage* image);
	void depthToPointCloud();
	void readDepthFiles(string &filename);
	void savePointCloud(int iteration, IplImage* image);
	void printCurrentGlobal_R_T();
};
