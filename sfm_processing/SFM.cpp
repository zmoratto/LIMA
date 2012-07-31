#include "SFM.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void printDepthWarning();

SFM::SFM(char* configFile, char* inputFilename)
{
	prevImage = NULL;
	prevDepth = NULL;
	currDepth = NULL;
	sfmInit = 1;
	string configFileString = string(configFile);
	string inputFilenameString = string(inputFilename);

	//Read Configuration Files and Set Up Tiles
	readConfigurationFile(configFileString);
	readImageFilenames(inputFilenameString);

	IplImage* image = cvLoadImage(inputFiles[configParams.firstFrame].c_str(), CV_LOAD_IMAGE_UNCHANGED);

	imageWidth = image->width;
	imageHeight = image->height;


	referenceTile.setUpReferenceTiles(image);

	numTiles = referenceTile.tiles.size();

	//Read Camera/Point Information Depending on Depth Info
	if (configParams.depthInfo == KINECT_DEPTH && depthFiles.size() == 0)
	{
		int cameraCalibrationFileReadError= readCameraCalibrationFile();
		pose.readDepthFiles(configParams.kinectDepthFilename, depthFiles);
	}
	if (configParams.depthInfo == STEREO_DISP && sfmInit)
	{
		int stereoCalibrationFileReadError = readStereoCalibrationFile();
	}
	if(configParams.depthInfo == POINT_CLOUD && depthFiles.size() == 0)
	{
		pose.readDepthFiles(configParams.pointCloudFilename, depthFiles);
	}

	//Set up 3D Point Vectors
	pose.resizeXYZVectors(imageWidth,imageHeight);

	cvReleaseImage(&image);
}

SFM::~SFM()
{
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
		fin >> identifier >> configParams.depthInfo;
		fin >> identifier >> configParams.showResultsFlag;
		fin >> identifier >> configParams.saveResultsFlag;
		fin >> identifier >> configParams.cameraCalibrationFilename;
		fin >> identifier >> configParams.stereoCalibrationFilename;
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
		fin.close();
	}
	else 
	{
	    printWarning(configurationFilename);
		restoreDefaultParameters();
	}  
}

void SFM::printConfigParams()
{
	cout<<"firstFrame="<<configParams.firstFrame<<endl;
	cout<<"lastFrame="<<configParams.lastFrame<<endl;
	cout<<"depthInfo="<<configParams.depthInfo<<endl;
	cout<<"showResultsFlag="<<configParams.showResultsFlag<<endl;
	cout<<"saveResultsFlag="<<configParams.saveResultsFlag<<endl;
	cout<<"cameraCalibrationFilename="<<configParams.cameraCalibrationFilename<<endl;
	cout<<"stereoCalibrationFilename="<<configParams.stereoCalibrationFilename<<endl;
	cout<<"minMatches="<<pose.minMatches<<endl;
	cout<<"detThresh="<<pose.detThresh<<endl;
	cout<<"featureMethod="<<feat.featureMethod<<endl;
	cout<<"poseEstimationMethod="<<configParams.poseEstimationType<<endl;
	cout<<"nndrRatio="<<pose.nndrRatio<<endl;
	cout<<"numNN="<<pose.numNN<<endl;
	cout<<"zDist="<<pose.zDist<<endl;
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
}


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
		camera_matrix = (cv::Mat_<double>(3,3)<<a00, a01, a02, a10, a11, a12, a20, a21, a22);

		stringstream sline1; 
		getline (calibrationFile,line);
		sline1 << line;
		sline1 >> identifier >> a00 >> a01 >> a02 >> a10;
		dist_coeffs = (cv::Mat_<double>(4,1)<<a00, a01, a02, a10);
	 
		calibrationFile.close();
		return 1;
	}
	else
	{
		cout << "Unable to open file " << configParams.cameraCalibrationFilename << endl; 
		return 0;
	}
}

int SFM::readStereoCalibrationFile()
{
}

void SFM::restoreDefaultParameters()
{
	configParams.firstFrame = 100;
	configParams.lastFrame = 1200;
	configParams.depthInfo = 0;
	configParams.showResultsFlag = 1;
	configParams.saveResultsFlag = 1;
	configParams.cameraCalibrationFilename = "camera_calibration.txt";
	configParams.stereoCalibrationFilename = "stereo_calibration.txt";
	configParams.pointCloudFilename = "rectified_left_list.txt";
	pose.minMatches = 8;
	pose.detThresh = 0.0000001;
	feat.featureMethod = 0;
	configParams.poseEstimationType = 0;
	pose.nndrRatio = 0.5;
	pose.numNN = 2;
	pose.zDist = 100;
	pose.xDist = 1000;
	pose.yDist = 1000;
	pose.nbMatches = 8;
	pose.homographyMethod = 8; //RANSAC
	pose.flannCheck = 1024;
	pose.ransacPix = 50;
	pose.ransacAccuracy = 0.99;
}

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

void SFM::processTile(IplImage *image, int frameIndex)
{
	int index, cnt = 0;
	pose.clearMatchedPixels();
	vector<KeyPoint> temp;
	temp.clear();
	vector<int> loc;
	Mat tempMat;
	int j = 0;	
	int left = imageWidth, right = 0, top = imageHeight, bottom = 0, x, y;
	CvRect mask;
	IplImage* imageMask;
	int tileY, tileX;
	KeyPoint tilePt;
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
	clock_t t1, t2;
	t1 = clock();
	feat.process(copy);
	t2 = clock();
	cout << "Feature Extraction Time: " << (double(t2)-double(t1))/CLOCKS_PER_SEC << endl;

	cvResetImageROI(copy);

	for(int i=0; i<feat.key_points.size(); i++)
	{
		feat.key_points.at(i).pt.x = feat.key_points.at(i).pt.x + double(mask.x);
		feat.key_points.at(i).pt.y = feat.key_points.at(i).pt.y + double(mask.y);
	}

	if(pose.globalCurrDescriptors.cols == 0)
		pose.globalCurrDescriptors.create(0, feat.point_descriptors.cols, feat.point_descriptors.type());

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
	tempMat.create(temp.size(), feat.point_descriptors.cols, feat.point_descriptors.type());

	cnt = 0;
	for(int i=0; i<feat.key_points.size(); i++)
	{
		index = (imageWidth)*int(feat.key_points.at(i).pt.y + tileY) + int(feat.key_points.at(i).pt.x + tileX);
		if(pose.curr_z[index] > 0)
		{
			Mat dest(tempMat.rowRange(cnt, cnt+1));
			feat.point_descriptors.row(i).copyTo(dest);
			cnt++;
		}
	}

	feat.point_descriptors.release();
	feat.key_points.clear();

	feat.point_descriptors = tempMat.clone();
	feat.key_points = temp;

	//FLANN Formatting
	if(feat.point_descriptors.type()!=CV_32F) 
		feat.point_descriptors.convertTo(pose.currObjectData, CV_32F);
	else 
		pose.currObjectData = feat.point_descriptors;

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
	Mat dest1(pose.globalCurrDescriptors.rowRange(oldRow, oldRow+feat.key_points.size()));
	pose.currObjectData.copyTo(dest1);
}

void SFM::preprocessing(IplImage* image, int frameIndex)
{
	int width = image->width;
	int height = image->height;

	//Read and Save Point Cloud File
	if (configParams.depthInfo == KINECT_DEPTH)
	{
		currDepth = cvLoadImage(depthFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);
		pose.depthToPointCloud(currDepth, camera_matrix);
	}
	else if (configParams.depthInfo == STEREO_DISP)
	{
		currDepth = cvLoadImage(depthFiles[frameIndex].c_str(), CV_LOAD_IMAGE_UNCHANGED);
		Mat currDispMat = &(*currDepth);
		pose.dispToPointCloud(currDispMat, Q);
		currDispMat.release();
	}
	else if(configParams.depthInfo == POINT_CLOUD)
	{
		pose.savePointCloud(depthFiles, frameIndex, image);
	}
}

void SFM::process(IplImage* image, int frameIndex)
{
	clock_t t1, t2;
	t1 = clock();
	preprocessing(image, frameIndex); //Read PC
	t2 = clock();
	referenceTile.process(image);

	cout << "PreProcessing Time: " << (double(t2)-double(t1))/CLOCKS_PER_SEC << endl;

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
		showGlobalMatches(prevImage, image, frameIndex);
	}

	pose.savePointProj(frameIndex, configParams.firstFrame);


	//Update
	update(image);

	t2 = clock();

	cout << "Overall Time: " << (double(t2)-double(t1))/CLOCKS_PER_SEC << endl;
}


void SFM::update(IplImage* image)
{
	pose.prevMatchedPoints.clear();
	pose.currMatchedPoints.clear();
	pose.prevMatchedPixels.clear();
	pose.currMatchedPixels.clear();
	pose.prevObjectData.release();
	pose.currObjectData.release();
	pose.prevKeypoints.clear();

	//Global Parameters
	pose.globalPrevMatchedPts.clear();
	pose.globalCurrMatchedPts.clear();
	pose.globalPrevMatchedPix.clear();
	pose.globalCurrMatchedPix.clear();
	pose.globalPrevKeyPoints = pose.globalCurrKeyPoints;
	pose.globalCurrKeyPoints.clear();
	pose.globalPrevDescriptors = pose.globalCurrDescriptors;
	pose.globalCurrDescriptors.release();

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

	cvReleaseImage(&image);
	image = NULL;
}

void SFM::clear()
{
	sfmInit = 1;
	feat.clear();
	pose.clear();
	if(prevImage != NULL)
		cvReleaseImage(&prevImage);
	prevImage = NULL;
}

void SFM::showGlobalMatches(IplImage* image1, IplImage* image2, int frameIndex)
{
	int index;
	IplImage* composedImage;
	int validIndex = 0;
	int validSize = pose.validPix.size();

	//Show all keypoints
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
	//Show KeyPoints Image2
	for(int m=0; m<pose.globalCurrKeyPoints.size(); m++)
	{
		index = m;
		pt1.x = pose.globalCurrKeyPoints[index].pt.x-1;
		pt2.x = pose.globalCurrKeyPoints[index].pt.x+1; 
		pt1.y = pose.globalCurrKeyPoints[index].pt.y-1;
		pt2.y = pose.globalCurrKeyPoints[index].pt.y+1;
		cvRectangle(image2, pt1, pt2, CV_RGB(0, 0, 0), 3, 8, 0 );
	}


	//Show KeyPoints Image 1
	//CvPoint pt1, pt2;
	for(int m=0; m<pose.globalPrevMatchedPix.size(); m++)
	{
		index = m;
		pt1.x = pose.globalPrevMatchedPix[index].x-1;
		pt2.x = pose.globalPrevMatchedPix[index].x+1; 
		pt1.y = pose.globalPrevMatchedPix[index].y-1;
		pt2.y = pose.globalPrevMatchedPix[index].y+1;
		cvRectangle(image1, pt1, pt2, CV_RGB(255, 255, 255), 3, 8, 0 );
	}
	//Show KeyPoints Image2
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


void SFM::cleanUp(IplImage*& im1, IplImage*& im2)
{
	pose.globalCurrMatchedPts.clear();
	pose.globalCurrMatchedPix.clear();
	pose.globalCurrKeyPoints.clear();
	pose.globalPrevMatchedPts.clear();
	pose.globalPrevMatchedPix.clear();
	pose.globalPrevKeyPoints.clear();
	pose.globalPrevDescriptors.release();
	pose.globalCurrDescriptors.release();
	cvReleaseImage(&im1);
	cvReleaseImage(&im2);
	cvReleaseImage(&prevImage);
	im1 = NULL;
	im2 = NULL;
	prevImage = NULL;
}

