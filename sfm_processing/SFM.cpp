#include "SFM.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void printDepthWarning();

SFM::SFM(char* configFile)
{
	prevImage = NULL;
	currDepth = NULL;
	sfmInit = 1;
	string configFileString = string(configFile);

	//Read Configuration Files
	readConfigurationFile(configFileString);
}

SFM::~SFM()
{
/*	cameraMatrix.release();
	distCoeffs.release();

	if(prevImage != NULL)
	{
		cvReleaseImage(&prevImage);
		prevImage = NULL;
	}
	if(currDepth != NULL)
	{
		cvReleaseImage(&currDepth);
		currDepth = NULL;
	}

	pointFiles.clear();

	for(int i=0; i<numTiles; i++)
	{
		cvReleaseImage(&referenceTile.tileImages[i]);
		referenceTile.tileImages[i] = NULL;
	}	

	inputFiles.clear();
	depthFiles.clear();
*/
}

void SFM::setUpSFM(char* inputFilename, IplImage* image)
{
	if(configParams.depthInfo != STEREO_DEPTH)
	{
		string inputFilenameString = string(inputFilename);
		readImageFilenames(inputFilenameString);
		image = cvLoadImage(inputFiles[configParams.firstFrame].c_str(), CV_LOAD_IMAGE_UNCHANGED);
	}

	imageWidth = image->width;
	imageHeight = image->height;

	//Set up tiles
	referenceTile.setUpReferenceTiles(image);
	numTiles = referenceTile.tiles.size();

	//Read Camera/Point Information Depending on Depth Info
	if (configParams.depthInfo == KINECT_DEPTH && depthFiles.size() == 0)
	{
		int cameraCalibrationFileReadError= readCameraCalibrationFile();
		readDepthFiles(configParams.kinectDepthFilename);
	}
	if(configParams.depthInfo == POINT_CLOUD && depthFiles.size() == 0)
	{
		readDepthFiles(configParams.pointCloudFilename);
	}

	//Set up 3D Point Vectors
	pose.allocatePCMemory(imageWidth,imageHeight);

	cvReleaseImage(&image);
	image = NULL;
}

void SFM::printUsage()
{
	cout << endl;
	cout << "******************** USAGE ********************" << endl;
	cout << "./sfm_test <configFile> <inputFile> <outputDir>" << endl;
	cout << "***********************************************" << endl;
	cout << endl;
}

void SFM::printWarning(string& filename)
{
	cout << endl;
	cout << "*****************************************************************" << endl;
	cout << "WARNING: Unable to open file "<<filename<<"..."<<endl; 
	cout << "Using default parameters..." << endl;
	cout << "*****************************************************************" << endl;
	cout << endl;
}

void SFM::readConfigurationFile(string& configurationFilename)
{
	ifstream fin (configurationFilename.c_str());
	string line;

	if(fin.is_open())
	{     
		string identifier;

		fin >> identifier >> configParams.firstFrame;
		fin >> identifier >> configParams.lastFrame;
		fin >> identifier >> configParams.frameStep;
		fin >> identifier >> configParams.depthInfo;
		fin >> identifier >> configParams.showResultsFlag;
		fin >> identifier >> configParams.saveResultsFlag;
		fin >> identifier >> configParams.cameraCalibrationFilename;
		fin >> identifier >> configParams.pointCloudFilename;
		fin >> identifier >> configParams.kinectDepthFilename;
		fin >> identifier >> pose.minMatches;
		fin >> identifier >> pose.detThresh;
		fin >> identifier >> feat.featureMethod;
		fin >> identifier >> configParams.poseEstimationType;
		fin >> identifier >> pose.nndrRatio;
		fin >> identifier >> pose.numNN;
		fin >> identifier >> pose.zDist;
		fin >> identifier >> pose.xDist;
		fin >> identifier >> pose.yDist;
		fin >> identifier >> pose.nbMatches;
		fin >> identifier >> pose.homographyMethod;
		fin >> identifier >> pose.flannCheck;
		fin >> identifier >> pose.ransacPix;
		fin >> identifier >> pose.ransacAccuracy;
		fin >> identifier >> referenceTile.tileParams.tileWidth;
		fin >> identifier >> referenceTile.tileParams.tileHeight;
		fin >> identifier >> referenceTile.tileParams.xOverlap;
		fin >> identifier >> referenceTile.tileParams.yOverlap;
		fin >> identifier >> referenceTile.tileParams.tileScaleX;
		fin >> identifier >> referenceTile.tileParams.tileScaleY;
		fin >> identifier >> configParams.pointProjFilename;
		fin >> identifier >> pose.weightedPortion;
		fin >> identifier >> configParams.resultImageFilename;
		fin >> identifier >> feat.hessianThresh;
		fin >> identifier >> feat.nOctaves;
		fin >> identifier >> feat.nOctaveLayers;
		fin >> identifier >> feat.extended;
		fin >> identifier >> feat.upright;
		fin >> identifier >> feat.octaves;
		fin >> identifier >> feat.octaveLayers;
		fin >> identifier >> pose.matchingMethod;
		fin.close();
	}
	else //If the config file doesn't open, resort to default parameters 
	{
	    printWarning(configurationFilename);
		restoreDefaultParameters();
	}  
}

void SFM::printConfigParams()
{
	cout<<"firstFrame="<<configParams.firstFrame<<endl;
	cout<<"lastFrame="<<configParams.lastFrame<<endl;
	cout<<"frameStep="<<configParams.frameStep<<endl;
	cout<<"depthInfo="<<configParams.depthInfo<<endl;
	cout<<"showResultsFlag="<<configParams.showResultsFlag<<endl;
	cout<<"saveResultsFlag="<<configParams.saveResultsFlag<<endl;
	cout<<"cameraCalibrationFilename="<<configParams.cameraCalibrationFilename<<endl;
	cout<<"pointCloudFilename="<<configParams.pointCloudFilename<<endl;
	cout<<"kinectDepthFilename="<<configParams.kinectDepthFilename;
	cout<<"minMatches="<<pose.minMatches<<endl;
	cout<<"detThresh="<<pose.detThresh<<endl;
	cout<<"featureMethod="<<feat.featureMethod<<endl;
	cout<<"poseEstimationMethod="<<configParams.poseEstimationType<<endl;
	cout<<"nndrRatio="<<pose.nndrRatio<<endl;
	cout<<"numNN="<<pose.numNN<<endl;
	cout<<"zDist="<<pose.zDist<<endl;
	cout<<"xDist="<<pose.xDist<<endl;
	cout<<"yDist="<<pose.yDist<<endl;
	cout<<"nbMatches="<<pose.nbMatches<<endl;
	cout<<"homographyMethod="<<pose.homographyMethod<<endl;
	cout<<"flannCheck="<<pose.flannCheck << endl;
	cout << "ransacPix=" << pose.ransacPix << endl;
	cout << "ransacAccuracy=" << pose.ransacAccuracy << endl;
	cout << "tileWidth=" << referenceTile.tileParams.tileWidth << endl;
	cout << "tileHeight=" << referenceTile.tileParams.tileHeight << endl;
	cout << "xOverlap=" << referenceTile.tileParams.xOverlap << endl;
	cout << "yOverlap=" << referenceTile.tileParams.yOverlap << endl;
	cout << "tileScaleX=" << referenceTile.tileParams.tileScaleX << endl;
	cout << "tileScaleY=" << referenceTile.tileParams.tileScaleY << endl;
	cout << "pointProjFilename=" << configParams.pointProjFilename << endl;
	cout << "weightedPortion=" << pose.weightedPortion << endl;
	cout << "resultImageFilename=" << configParams.resultImageFilename << endl;
	cout << "hessianThresh=" << feat.hessianThresh <<endl;
	cout << "nOctaves=" << feat.nOctaves << endl;
	cout << "nOctaveLayers=" << feat.nOctaveLayers << endl;
	cout << "extended=" << feat.extended << endl;
	cout << "upright=" << feat.upright << endl;
	cout << "octaves=" << feat.octaves << endl;
	cout << "octaveLayers=" << feat.octaveLayers << endl;
	cout << "matchingMethod=" << pose.matchingMethod << endl;
}

//Read camera calibration file for kinect
int SFM::readCameraCalibrationFile()
{
	ifstream calibrationFile (configParams.cameraCalibrationFilename.c_str());
	string line;
	double a00, a01,a02, a10, a11, a12, a20, a21, a22;

	if (calibrationFile.is_open())
	{     
		string identifier;

		stringstream sline;  
		getline (calibrationFile,line);
		sline<<line;
		sline >> identifier >> a00 >> a01 >> a02 >> a10 >> a11>> a12 >> a20 >> a21 >> a22;
		cameraMatrix = (cv::Mat_<double>(3,3)<<a00, a01, a02, a10, a11, a12, a20, a21, a22);

		stringstream sline1; 
		getline (calibrationFile,line);
		sline1 << line;
		sline1 >> identifier >> a00 >> a01 >> a02 >> a10;
		distCoeffs = (cv::Mat_<double>(4,1)<<a00, a01, a02, a10);
	 
		calibrationFile.close();
		return 1;
	}
	else
	{
		cout << "Unable to open file " << configParams.cameraCalibrationFilename << endl; 
		return 0;
	}
}

//Restore default parameters if no config file
void SFM::restoreDefaultParameters()
{
	configParams.firstFrame = 0;
	configParams.lastFrame = 10;
	configParams.frameStep = 1;
	configParams.depthInfo = 0;
	configParams.showResultsFlag = 1;
	configParams.saveResultsFlag = 1;
	configParams.cameraCalibrationFilename = "camera_calibration.txt";
	configParams.pointCloudFilename = "point_cloud_list.txt";
	configParams.kinectDepthFilename = "kinect_depth_list.txt";
	pose.minMatches = 8;
	pose.detThresh = 0.0000001;
	feat.featureMethod = 0;
	configParams.poseEstimationType = 0;
	pose.nndrRatio = 0.55;
	pose.numNN = 2;
	pose.zDist = 800;
	pose.xDist = 1000;
	pose.yDist = 1000;
	pose.nbMatches = 8;
	pose.homographyMethod = 8;
	pose.flannCheck = 64;
	pose.ransacPix = 3;
	pose.ransacAccuracy = 0.99;
	referenceTile.tileParams.tileWidth = 500;
	referenceTile.tileParams.tileHeight = 500;
	referenceTile.tileParams.xOverlap = 0;
	referenceTile.tileParams.yOverlap = 0;
	referenceTile.tileParams.tileScaleX = 0;
	referenceTile.tileParams.tileScaleY = 0;
	configParams.pointProjFilename = "results/pointProj.txt";
	pose.weightedPortion = 0.5;
	configParams.resultImageFilename = "results/composedImage";
	feat.hessianThresh = 800.0;
	feat.nOctaves = 4;
	feat.nOctaveLayers = 2;
	feat.extended = 0;
	feat.upright = 1;
	feat.octaves = 3;
	feat.octaveLayers = 4;
	pose.matchingMethod = 1;
}

//Read list of input images and save into vector
void SFM::readImageFilenames(string& inputFile)
{
	string image, depth;
	int i = 0;
	ifstream fin;
	string line;

	fin.open(inputFile.c_str());

	if(!fin.is_open())
	{
		cout << "ERROR: Unable to open file " << inputFile << endl;
	}
	while(fin.good())
	{
		stringstream ss;
		getline(fin, line);
		ss << line;
		ss >> image;
		inputFiles.push_back(image);
	}
	fin.close();
}

//Perform feature extraction on tile and save extracted features and descriptors in global collection variables
void SFM::processTile(IplImage *image, int frameIndex)
{
	int index, cnt = 0;
	pose.currMatchedPixels.clear();
	pose.prevMatchedPixels.clear();
	vector<cv::KeyPoint> temp;
	temp.clear();
	vector<int> loc;
	cv::Mat tempMat;
	int j = 0;	
	int left = imageWidth, right = 0, top = imageHeight, bottom = 0, x, y;
	CvRect mask;
	IplImage* imageMask;
	int tileY, tileX;
	cv::KeyPoint tilePt;
	int oldRow = 0;
	int buffer = 25;
	IplImage* tempImg = cvCreateImage(cvSize(imageWidth,imageHeight),image->depth,image->nChannels);
	IplImage* copy = cvCreateImage(cvGetSize(image),image->depth,image->nChannels);

	//Clear Everything
	feat.key_points.clear();
	feat.point_descriptors.release();
	pose.prevKeypoints.clear();
	pose.currMatchedPoints.clear();
	pose.prevMatchedPoints.clear();
	pose.prevMatchedPixels.clear();
	pose.currMatchedPixels.clear();

	//Set up mask for FE
	//Only extract features from area where 3D information is available
	for(int i=0; i<imageWidth*imageHeight; i++)
	{
		if(pose.curr_z[i] > 0)
		{
			y = i/imageWidth;
			x = i%imageWidth;
			if(x < left)
				left = x;
			if(x > right)
				right = x;
			if(y < top)
				top = y;
			if(y > bottom)
				bottom = y;
		}
	}

	mask.x = referenceTile.tiles[currTile].x;
	mask.y = referenceTile.tiles[currTile].y;
	mask.width = referenceTile.tiles[currTile].width;
	mask.height = referenceTile.tiles[currTile].height;

	if(mask.x < left)
		mask.x = left;
	if(referenceTile.tiles[currTile].x + referenceTile.tiles[currTile].width > right)	
		mask.width = right - mask.x;

	if(mask.y < top)
		mask.y = top;
	if(referenceTile.tiles[currTile].y + referenceTile.tiles[currTile].height > bottom)
		mask.height = bottom - mask.y;

	mask.x -= referenceTile.tiles[currTile].x;
	mask.y -= referenceTile.tiles[currTile].y;


	cvCopy(image, copy, NULL);
	cvSetImageROI(copy, mask);

	//Feature Extraction
	feat.process(copy);
	cvResetImageROI(copy);

	//FLANN Formatting
	if(feat.point_descriptors.type()!=CV_32F) 
		feat.point_descriptors.convertTo(pose.currObjectData, CV_32F);
	else 
		pose.currObjectData = feat.point_descriptors;

	for(int i=0; i<feat.key_points.size(); i++)
	{
		feat.key_points.at(i).pt.x = feat.key_points.at(i).pt.x + double(mask.x);
		feat.key_points.at(i).pt.y = feat.key_points.at(i).pt.y + double(mask.y);
	}

	if(pose.globalCurrDescriptors.cols == 0)
		pose.globalCurrDescriptors.create(0, pose.currObjectData.cols, CV_32F);

	//Remove KeyPoints with No Depth
	tileY = referenceTile.tiles[currTile].y;
	tileX = referenceTile.tiles[currTile].x;
	for(int i=0; i<feat.key_points.size(); i++)
	{
		index = (imageWidth)*int(feat.key_points.at(i).pt.y + tileY) + int(feat.key_points.at(i).pt.x + tileX);
		if(pose.curr_z[index] > 0)
		{
			temp.push_back(feat.key_points.at(i));
			loc.push_back(i);
		}
	}
	tempMat.create(temp.size(), pose.currObjectData.cols, CV_32F);

	cnt = 0;
	for(int i=0; i<feat.key_points.size(); i++)
	{
		index = (imageWidth)*int(feat.key_points.at(i).pt.y + tileY) + int(feat.key_points.at(i).pt.x + tileX);
		if(pose.curr_z[index] > 0)
		{
			cv::Mat dest(tempMat.rowRange(cnt, cnt+1));
			pose.currObjectData.row(i).copyTo(dest);
			cnt++;
		}
	}

	pose.currObjectData.release();
	feat.key_points.clear();

	pose.currObjectData = tempMat.clone();
	feat.key_points = temp;

	for(int i=0; i<feat.key_points.size(); i++)
	{
		tilePt = feat.key_points[i];
		tilePt.pt.x = feat.key_points[i].pt.x + referenceTile.tiles[currTile].x;
		tilePt.pt.y = feat.key_points[i].pt.y + referenceTile.tiles[currTile].y;
		pose.globalCurrKeyPoints.push_back(tilePt);
	}

	//Collect All Descriptors
	oldRow = pose.globalCurrDescriptors.rows;
	pose.globalCurrDescriptors.resize(oldRow+feat.key_points.size());
	cv::Mat dest1(pose.globalCurrDescriptors.rowRange(oldRow, oldRow+feat.key_points.size()));
	pose.currObjectData.copyTo(dest1);

	cvReleaseImage(&copy);
}

void SFM::preprocessing(IplImage* image, int frameIndex)
{
	int width = image->width;
	int height = image->height;

	//Read and Save Point Cloud File
	if (configParams.depthInfo == KINECT_DEPTH)
	{
		currDepth = cvLoadImage(depthFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);
		depthToPointCloud();
	}
	else if(configParams.depthInfo == POINT_CLOUD)
	{
		savePointCloud(frameIndex, image);
	}
}

//Process each tile, process pose and map point clouds
void SFM::process(IplImage* image, int frameIndex)
{
	clock_t t1, t2;
	preprocessing(image, frameIndex); //Read PC
	referenceTile.process(image);

	t1 = clock();

	//Process Each Tile
	for(int j=0; j<numTiles;j++)
	{
		currTile = j;
		processTile(referenceTile.tileImages[j], frameIndex);
	}

	pose.removeDuplicates(image);

	if(configParams.firstFrame != frameIndex)
	{
		pose.process(configParams.depthInfo, image);
		t2 = clock();
		cout << "Overall Time: " << (double(t2)-double(t1))/CLOCKS_PER_SEC << endl;

		if(configParams.saveResultsFlag)
		{
			printCurrentGlobal_R_T();
			string resultsString = string("results/PointCloud.txt");
			mapping(resultsString, image);
			showGlobalMatches(prevImage, image, frameIndex);
		}
	}

	pose.savePointProj(frameIndex, configParams.firstFrame);

	//Update
	update(image);
}

//Move current points and matches to previous variables
void SFM::update(IplImage* image)
{
	pose.clear();

	//Global Parameters
	pose.globalPrevMatchedPts.clear();
	pose.globalCurrMatchedPts.clear();
	pose.globalPrevMatchedPix.clear();
	pose.globalCurrMatchedPix.clear();
	pose.globalPrevKeyPoints = pose.globalCurrKeyPoints;
	pose.globalCurrKeyPoints.clear();
	pose.globalPrevDescriptors = pose.globalCurrDescriptors.clone();
	pose.globalCurrDescriptors.release();
	pose.globalCurrDescriptors.create(0, pose.globalCurrDescriptors.cols, CV_32F);

	//3D Points
	pose.resetPrevPos();

	if(prevImage != NULL)
	{
		cvReleaseImage(&prevImage);
		prevImage = NULL;
	}

	prevImage = cvCreateImage(cvGetSize(image), image->depth, image->nChannels);

	cvCopy(image, prevImage);

	if(currDepth != NULL)
	{
		cvReleaseImage(&currDepth);
		currDepth = NULL;
	}
	sfmInit = 0;

	referenceTile.clear();
	pose.matchWeights.clear();

	for(int i=0; i<numTiles; i++)
	{
		cvReleaseImage(&referenceTile.tileImages[i]);
		referenceTile.tileImages[i] = NULL;
	}

/*	if(image != NULL)
	{
		cvReleaseImage(&image);
		image = NULL;
	}
*/
}

//Creates composed image and displays matches and keypoints
void SFM::showGlobalMatches(IplImage* image1, IplImage* image2, int frameIndex)
{
	int index;
	IplImage* composedImage;
	int validIndex = 0;
	int validSize = pose.validPix.size();

	//Show All Key Points
	CvPoint pt1, pt2;
	for(int m=0; m<pose.globalPrevKeyPoints.size(); m++)
	{
		index = m;
		pt1.x = pose.globalPrevKeyPoints[index].pt.x-1;
		pt2.x = pose.globalPrevKeyPoints[index].pt.x+1; 
		pt1.y = pose.globalPrevKeyPoints[index].pt.y-1;
		pt2.y = pose.globalPrevKeyPoints[index].pt.y+1;
		cvRectangle(image1, pt1, pt2, CV_RGB(0, 0, 0), 3, 8, 0 );
	}
	for(int m=0; m<pose.globalCurrKeyPoints.size(); m++)
	{
		index = m;
		pt1.x = pose.globalCurrKeyPoints[index].pt.x-1;
		pt2.x = pose.globalCurrKeyPoints[index].pt.x+1; 
		pt1.y = pose.globalCurrKeyPoints[index].pt.y-1;
		pt2.y = pose.globalCurrKeyPoints[index].pt.y+1;
		cvRectangle(image2, pt1, pt2, CV_RGB(0, 0, 0), 3, 8, 0 );
	}


	//Matching Key Points
	for(int m=0; m<pose.globalPrevMatchedPix.size(); m++)
	{
		index = m;
		pt1.x = pose.globalPrevMatchedPix[index].x-1;
		pt2.x = pose.globalPrevMatchedPix[index].x+1; 
		pt1.y = pose.globalPrevMatchedPix[index].y-1;
		pt2.y = pose.globalPrevMatchedPix[index].y+1;
		cvRectangle(image1, pt1, pt2, CV_RGB(255, 255, 255), 3, 8, 0 );
	}
	for(int m=0; m<pose.globalPrevMatchedPix.size(); m++)
	{
		index = m;
		pt1.x = pose.globalCurrMatchedPix[index].x-1;
		pt2.x = pose.globalCurrMatchedPix[index].x+1; 
		pt1.y = pose.globalCurrMatchedPix[index].y-1;
		pt2.y = pose.globalCurrMatchedPix[index].y+1;
		cvRectangle(image2, pt1, pt2, CV_RGB(255, 255, 255), 3, 8, 0 );
	}

	composedImage = cvCreateImage(cvSize(image1->width*2, image1->height), image1->depth, image1->nChannels);
	CvRect imageRect = cvRect(0,0, image1->width, image1->height); 
	cvSetImageROI(composedImage, imageRect);
	cvCopy(image2, composedImage);
	cvResetImageROI(composedImage);
						  
	imageRect = cvRect((image2->width)-1, 0, image2->width, image2->height);
	cvSetImageROI(composedImage, imageRect);
	cvCopy(image1, composedImage);
	cvResetImageROI(composedImage);

	//Show Image matches
	for(int m=0; m<pose.globalPrevMatchedPix.size(); m++)
	{
		index = m;
		pt1.x = pose.globalCurrMatchedPix[index].x;
		pt1.y = pose.globalCurrMatchedPix[index].y;
		if ((pt1.x <= 0) || (pt1.y <=0))
		{
			cout<<"wrong: x="<<pt1.x<<", y="<<pt1.y << "  Index=" << index << "  m=" << m <<endl;
		}
		pt2.x = pose.globalPrevMatchedPix[index].x+(composedImage->width/2);
		pt2.y = pose.globalPrevMatchedPix[index].y;

		if(configParams.depthInfo != 0)
		{
			if(pose.validPix.size() != 0 && pose.validPix[validIndex] == index)
			{
				cvLine( composedImage, pt1, pt2, CV_RGB(255,255,255), 1);
				if(validSize != validIndex-1)
					validIndex++;
			}
			else
				cvLine( composedImage, pt1, pt2, CV_RGB(0,0,0), 1);
		}
		else
		{
			cvLine( composedImage, pt1, pt2, CV_RGB(255,255,255), 1);
		}
	}

	stringstream ss;
	ss<<frameIndex;
	string composedFilename = configParams.resultImageFilename + ss.str() + string(".jpg");
	cvSaveImage(composedFilename.c_str(), composedImage);
	cvReleaseImage(&composedImage);
}

//Takes global R and T and applies them to the point clouds to get aligned point clouds.
void SFM::mapping(string& filename, IplImage* image)
{
	ofstream fout;
	stringstream ss;
	int width = image->width;
	float rgb;
	uchar* data = (uchar*)image->imageData;
	uchar *rData, *gData, *bData;
	CvMat* temp_T = cvCreateMat(3,1,CV_32FC1);
	CvMat* temp_R = cvCreateMat(3,3,CV_32FC1);
	int numChannels = image->nChannels;

	//Split the image into 3 separate channels to get the colors for the point cloud visualization.
	IplImage* r = cvCreateImage( cvGetSize(image), image->depth,1 );
	IplImage* g = cvCreateImage( cvGetSize(image), image->depth,1 );
	IplImage* b = cvCreateImage( cvGetSize(image), image->depth,1 );

	if(numChannels == 3)
	{
		cvSplit(image,b,g,r,NULL);
		rData = (uchar*)r->imageData;
		gData = (uchar*)g->imageData;
		bData = (uchar*)b->imageData;
	}

	int numPoints = pose.curr_x.size();
	int counter = 0;
	CvMat *pointMat = cvCreateMat (3, 1, CV_32FC1); //Used to hold each 3D point for the matrix multiplication
	CvMat *outPointMat = cvCreateMat (3, 1, CV_32FC1); //Holds the output of each 3D point matrix multiplication

	//Compute the opposite of globalT to apply to the point clouds
	temp_T->data.fl[0] = -1*pose.currGlobal_T->data.fl[0];
	temp_T->data.fl[1] = -1*pose.currGlobal_T->data.fl[1];
	temp_T->data.fl[2] = -1*pose.currGlobal_T->data.fl[2];

	//Compute the transpose of globalR to apply to the point clouds
	temp_R->data.fl[0] = pose.currGlobal_R->data.fl[0];
	temp_R->data.fl[1] = pose.currGlobal_R->data.fl[3];
	temp_R->data.fl[2] = pose.currGlobal_R->data.fl[6];
	temp_R->data.fl[3] = pose.currGlobal_R->data.fl[1];
	temp_R->data.fl[4] = pose.currGlobal_R->data.fl[4];
	temp_R->data.fl[5] = pose.currGlobal_R->data.fl[7];
	temp_R->data.fl[6] = pose.currGlobal_R->data.fl[2];
	temp_R->data.fl[7] = pose.currGlobal_R->data.fl[5];
	temp_R->data.fl[8] = pose.currGlobal_R->data.fl[8];

	for (int i = 0; i < numPoints; i++)
	{
		if(pose.curr_z[i] > 0) //If the point is in front of the camera, apply the R and T to the point
		{
			//Save the current 3D point into pointMat
			cvSetReal2D(pointMat, 0, 0, pose.curr_x[i]);
			cvSetReal2D(pointMat, 1, 0, pose.curr_y[i]);
			cvSetReal2D(pointMat, 2, 0, pose.curr_z[i]);

			//Apply globalR Transpose and negative globalT to the current point
			cvMatMulAdd(temp_R, pointMat, temp_T, outPointMat);

			//Write out each mapped point and the color information associated with it.
			if(numChannels == 1)
				ss << outPointMat->data.fl[0] << " " << outPointMat->data.fl[1] << " " << outPointMat->data.fl[2] << " " << int(data[i]) << " " << int(data[i]) << " " << int(data[i]) << '\n';
			else if(numChannels == 3)
			{
				ss << outPointMat->data.fl[0] << " " << outPointMat->data.fl[1] << " " << outPointMat->data.fl[2] << " " << int(rData[i]) << " " << int(gData[i]) << " " << int(bData[i]) << '\n';
			}
			counter++;

		}
	}

	//Clean up
	cvReleaseMat(&pointMat);
	cvReleaseMat(&outPointMat);
	cvReleaseImage(&r);
	cvReleaseImage(&g);
	cvReleaseImage(&b);
	fout.open (filename.c_str(), ios::app);
	fout << ss.str();
	fout.close();
	ss.str("");
}

//Compute depth from kinect depth image
void SFM::depthToPointCloud()
{
	const int GAMMASIZE =  2048;
	float gamma[GAMMASIZE];
	const float k1 = 1.1863;
	const float k2 = 2842.5;
	const float k3 = 0.1236;

	for (size_t i = 0; i < GAMMASIZE; i++)
		gamma[i] = k3 * tan(i/k2 + k1);

	// camera intrinsic parameters, representative values, see http://nicolas.burrus.name/index.php/Research/KinectCalibration for more info
	float cx = cameraMatrix.at<double>(0,2);
	float cy = cameraMatrix.at<double>(1,2);
	float fx = cameraMatrix.at<double>(0,0);
	float fy = cameraMatrix.at<double>(1,1);

	unsigned char *depthData = (unsigned char*)(currDepth->imageData);
	for (int i = 0; i < currDepth->height; i++)
	{
		for (int j = 0; j < currDepth->width; j++)
		{
			float gamma_depth_data = gamma[((short*)depthData)[i*currDepth->width+j]];
			int index = i*currDepth->width+j;
			pose.curr_x[index] = (j - cx) * gamma_depth_data / fx;
			pose.curr_y[index] = (i - cy) * gamma_depth_data / fy;
			pose.curr_z[index] = gamma_depth_data;
		}
	}
}

//Reads depth filenames from input file (kinect and point clouds)
void SFM::readDepthFiles(string& filename)
{
	string temp;
	ifstream fin;

	fin.open(filename.c_str());

	while(fin.good())
	{
		fin >> temp;
		depthFiles.push_back(temp);
	}

	fin.close();
}

//Saves individual point clouds from the file into memory
void SFM::savePointCloud(int iteration, IplImage* image)
{
	ifstream fin;
	string filename;
	int index = 0;
	float x, y, z;
	int width = image->width;
	float xPos, yPos;
	stringstream ss;
	ss<<iteration;
	filename = depthFiles[iteration];
	fin.open(filename.c_str());
	string line;

	for(int i=0; i<image->width*image->height; i++)
	{
		pose.curr_x[i] = 0;
		pose.curr_y[i] = 0;
		pose.curr_z[i] = 0;
	}

	while(fin.good())
	{
		getline(fin, line);

		sscanf(line.c_str(), "%f %f %f %f %f", &yPos, &xPos, &x, &y, &z);

		index = yPos*width + xPos;

		pose.curr_x[index] = x;
		pose.curr_y[index] = y;
		pose.curr_z[index] = z;
	}

	fin.close();
}

void SFM::printCurrentGlobal_R_T()
{
	cout<<"global_R=["<<pose.currGlobal_R->data.fl[0]<<" "<<pose.currGlobal_R->data.fl[1]<<" "<<pose.currGlobal_R->data.fl[2]<<"]"<<endl;
	cout<<"         ["<<pose.currGlobal_R->data.fl[3]<<" "<<pose.currGlobal_R->data.fl[4]<<" "<<pose.currGlobal_R->data.fl[5]<<"]"<<endl;
	cout<<"         ["<<pose.currGlobal_R->data.fl[6]<<" "<<pose.currGlobal_R->data.fl[7]<<" "<<pose.currGlobal_R->data.fl[8]<<"]"<<endl;
	cout << endl;
	cout<<"global_T=["<<pose.currGlobal_T->data.fl[0]<<" "<<pose.currGlobal_T->data.fl[1]<<" "<<pose.currGlobal_T->data.fl[2]<<"]"<<endl;
}
