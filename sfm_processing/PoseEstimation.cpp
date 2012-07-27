#include <opencv2/core/core.hpp>
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include "pcl/common/common_headers.h"
#include <pcl/correspondence.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/registration/transformation_estimation_lm.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

#include <Eigen/Dense>

#include "PoseEstimation.h"

using namespace Eigen;
using namespace std;
using namespace pcl;

PoseEstimation::PoseEstimation()
{
	//Clear Point Vectors
	currMatchedPoints.clear();
	prevMatchedPoints.clear();

	//Previous Global
	prevGlobal_R = cvCreateMat(3,3,CV_32F);
	prevGlobal_T = cvCreateMat(3,1,CV_32F);

	//Initialize with first frame
	cvSetReal2D(prevGlobal_R, 0, 0, 1.0);
	cvSetReal2D(prevGlobal_R, 0, 1, 0.0);
	cvSetReal2D(prevGlobal_R, 0, 2, 0.0);
	cvSetReal2D(prevGlobal_R, 1, 0, 0.0);
	cvSetReal2D(prevGlobal_R, 1, 1, 1.0);
	cvSetReal2D(prevGlobal_R, 1, 2, 0.0);
	cvSetReal2D(prevGlobal_R, 2, 0, 0.0);
	cvSetReal2D(prevGlobal_R, 2, 1, 0.0);
	cvSetReal2D(prevGlobal_R, 2, 2, 1.0);

	cvSetReal2D(prevGlobal_T, 0, 0, 0.0);
	cvSetReal2D(prevGlobal_T, 1, 0, 0.0);
	cvSetReal2D(prevGlobal_T, 2, 0, 0.0);	

	//Current Global
	currGlobal_R = cvCreateMat(3,3,CV_32F);
	currGlobal_T = cvCreateMat(3,1,CV_32F);

	//Relative R & T
	relative_R = cvCreateMat (3, 3, CV_32FC1);
	relative_T = cvCreateMat(3,1,CV_32F);

	projectedPoints.clear();

}

PoseEstimation::~PoseEstimation()
{
	outlierMask.release();

	cvReleaseMat(&prevGlobal_R);
	cvReleaseMat(&prevGlobal_T);
	
	cvReleaseMat(&currGlobal_R);
	cvReleaseMat(&currGlobal_T);

	cvReleaseMat(&relative_R);
	cvReleaseMat(&relative_T);
}

int PoseEstimation::estimateRelativePose()
{
	cv::Point3f translation;
	cv::Point3f currCenter;
	cv::Point3f prevCenter;
	int i=0;

	int numValidMatches = 0;
	float numWeightedMatches = 0.0;
	int numMatches = globalCurrMatchedPts.size();

	validPix.clear();

	//Find Translation, Previous Center, and Current Center
	for (i = 0; i < numMatches; i++)
	{
		
		if ((outlierMask.at<bool>(i) == 1) && (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z) < zDist) && (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y) < yDist) && (fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x) < xDist))
		{
			translation = translation + matchWeights[i]*(globalCurrMatchedPts[i]-globalPrevMatchedPts[i]);
			if(matchWeights[i] != NO_WEIGHT)
			{
				currCenter = currCenter + (globalCurrMatchedPts[i]);
				prevCenter = prevCenter + (globalPrevMatchedPts[i]);
				validPix.push_back(i);
				numValidMatches++;
			}

			numWeightedMatches+=matchWeights[i];
		}
	}

	cout << "NumMatches: " << numMatches << endl;
	cout << "NumValidMatches: " << numValidMatches << endl;
	cout << "NumWeights: " << numWeightedMatches << endl;

	if (numValidMatches > 0)
	{
		translation.x = translation.x/numWeightedMatches;
		translation.y = translation.y/numWeightedMatches;
		translation.z = translation.z/numWeightedMatches;
		cvSetReal2D(relative_T, 0, 0, translation.x);
		cvSetReal2D(relative_T, 1, 0, translation.y);
		cvSetReal2D(relative_T, 2, 0, translation.z);

		currCenter.x = currCenter.x/numValidMatches;
		currCenter.y = currCenter.y/numValidMatches;
		currCenter.z = currCenter.z/numValidMatches;
		prevCenter.x = prevCenter.x/numValidMatches;
		prevCenter.y = prevCenter.y/numValidMatches;
		prevCenter.z = prevCenter.z/numValidMatches;

	}
	else
	{
		return 1; //Return ERROR
	}

	float a00 = 0.0;
	float a10 = 0.0; 
	float a20 = 0.0; 
	float a01 = 0.0; 
	float a11 = 0.0; 
	float a21 = 0.0; 
	float a02 = 0.0; 
	float a12 = 0.0; 
	float a22 = 0.0;

	//Compute Relative_R
	for (i = 0; i < numMatches; i++)
	{
		if ((outlierMask.at<bool>(i) == 1) && (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z) < zDist) && (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y) < yDist) && (fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x) < xDist))
		{
			if(matchWeights[i] != NO_WEIGHT)
			{
				a00 = a00 + (globalPrevMatchedPts[i].x-prevCenter.x)*(globalCurrMatchedPts[i].x - currCenter.x);
				a10 = a10 + (globalPrevMatchedPts[i].y-prevCenter.y)*(globalCurrMatchedPts[i].x - currCenter.x);
				a20 = a20 + (globalPrevMatchedPts[i].z-prevCenter.z)*(globalCurrMatchedPts[i].x - currCenter.x);

				a01 = a01 + (globalPrevMatchedPts[i].x-prevCenter.x)*(globalCurrMatchedPts[i].y - currCenter.y);
				a11 = a11 + (globalPrevMatchedPts[i].y-prevCenter.y)*(globalCurrMatchedPts[i].y - currCenter.y);
				a21 = a21 + (globalPrevMatchedPts[i].z-prevCenter.z)*(globalCurrMatchedPts[i].y - currCenter.y);

				a02 = a02 + (globalPrevMatchedPts[i].x-prevCenter.x)*(globalCurrMatchedPts[i].z - currCenter.z);
				a12 = a12 + (globalPrevMatchedPts[i].y-prevCenter.y)*(globalCurrMatchedPts[i].z - currCenter.z);
				a22 = a22 + (globalPrevMatchedPts[i].z-prevCenter.z)*(globalCurrMatchedPts[i].z - currCenter.z);
			}
		}
	}


	CvMat *AMat = cvCreateMat(3, 3, CV_32F);	
	cvSetReal2D(AMat, 0, 0, a00);
	cvSetReal2D(AMat, 1, 0, a10);
	cvSetReal2D(AMat, 2, 0, a20);
	cvSetReal2D(AMat, 0, 1, a01);
	cvSetReal2D(AMat, 1, 1, a11);
	cvSetReal2D(AMat, 2, 1, a21);
	cvSetReal2D(AMat, 0, 2, a02);
	cvSetReal2D(AMat, 1, 2, a12);
	cvSetReal2D(AMat, 2, 2, a22);

	if (numValidMatches >= minMatches)
	{

		CvMat *WMat = cvCreateMat (3, 1, CV_32FC1);
		CvMat *UMat = cvCreateMat (3, 3, CV_32FC1);
		CvMat *VMat = cvCreateMat (3, 3, CV_32FC1);

		cvSVD(AMat, WMat, UMat, VMat, CV_SVD_V_T);

		double det = cvDet(AMat);
		if (fabs(det)< detThresh)
		{
			cvReleaseMat(&VMat);
			cvReleaseMat(&WMat);
			cvReleaseMat(&UMat);
			cvReleaseMat(&AMat);
			return 1; //ERROR
		}
		else
		{
			if (det > 0)
			{
				cvGEMM(UMat, VMat, 1, NULL, 0, relative_R);
			}
			if (det < 0)
			{
				CvMat *signMat = cvCreateMat (3, 3, CV_32FC1);
				cvSetReal2D(signMat, 0, 0, 1);
				cvSetReal2D(signMat, 0, 1, 0);
				cvSetReal2D(signMat, 0, 2, 0);
	
				cvSetReal2D(signMat, 1, 0, 0);
				cvSetReal2D(signMat, 1, 1, 1);
				cvSetReal2D(signMat, 1, 2, 0);
	
				cvSetReal2D(signMat, 2, 0, 0);
				cvSetReal2D(signMat, 2, 1, 0);
				cvSetReal2D(signMat, 2, 2, -1);
				cvGEMM(UMat, signMat, 1, NULL, 0, relative_R);
				cvGEMM(relative_R, VMat, 1, NULL, 0, relative_R);

				cvReleaseMat(&signMat);	
			}
			cvReleaseMat(&VMat);
			cvReleaseMat(&WMat);
			cvReleaseMat(&UMat);
			cvReleaseMat(&AMat);
		}
		return 0;
	}
	else
	{
		cvReleaseMat(&AMat);
		return 1; //ERROR
	}
}


//ComposePose
void PoseEstimation::composePose()
{
	//Compute Global R and Global T using relative R and relative T.
	cvMatMul(prevGlobal_R, relative_R, currGlobal_R);
	cvMatMulAdd(prevGlobal_R, relative_T, prevGlobal_T, currGlobal_T);
}

void PoseEstimation::printCurrentGlobal_R_T()
{
	cout<<"global_R=["<<currGlobal_R->data.fl[0]<<" "<<currGlobal_R->data.fl[1]<<" "<<currGlobal_R->data.fl[2]<<"]"<<endl;
	cout<<"         ["<<currGlobal_R->data.fl[3]<<" "<<currGlobal_R->data.fl[4]<<" "<<currGlobal_R->data.fl[5]<<"]"<<endl;
	cout<<"         ["<<currGlobal_R->data.fl[6]<<" "<<currGlobal_R->data.fl[7]<<" "<<currGlobal_R->data.fl[8]<<"]"<<endl;
	cout << endl;
	cout<<"global_T=["<<currGlobal_T->data.fl[0]<<" "<<currGlobal_T->data.fl[1]<<" "<<currGlobal_T->data.fl[2]<<"]"<<endl;
}

void PoseEstimation::clearMatchedPoints()
{
	currMatchedPoints.clear();
	prevMatchedPoints.clear();
}

void PoseEstimation::clearMatchedPixels()
{
	currMatchedPixels.clear();
	prevMatchedPixels.clear();
}

void PoseEstimation::push_prev3Dpt(int index)
{
	cv::Point3f pt3;
	pt3.x = prev_x[index];
	pt3.y = prev_y[index];  
	pt3.z = prev_z[index];
	globalPrevMatchedPts.push_back(pt3);
}


void PoseEstimation::push_curr3Dpt(int index)
{
	cv::Point3f pt3;
	pt3.x = curr_x[index];
	pt3.y = curr_y[index];  
	pt3.z = curr_z[index];
	globalCurrMatchedPts.push_back(pt3);
}

int PoseEstimation::find_prev_index(int loc, int width)
{
	int prevCol = globalPrevMatchedPix[loc].x;
	int prevRow = globalPrevMatchedPix[loc].y;
	int prevIndex = prevRow*width + prevCol;

	return prevIndex;
}

int PoseEstimation::find_curr_index(int loc, int width)
{
	int currCol = globalCurrMatchedPix[loc].x;
	int currRow = globalCurrMatchedPix[loc].y;
	int currIndex = currRow*width + currCol;

	return currIndex;
}

void PoseEstimation::find_homography()
{
	cv::Mat H;

	if(homographyMethod == CV_RANSAC)
		H = cv::findFundamentalMat(globalCurrMatchedPix, globalPrevMatchedPix, outlierMask, CV_FM_RANSAC, ransacPix, ransacAccuracy);
	if(homographyMethod == CV_LMEDS)
		H = cv::findHomography(globalPrevMatchedPix, globalCurrMatchedPix, outlierMask, CV_LMEDS);
	if(homographyMethod == 0)
		H = cv::findHomography(globalPrevMatchedPix, globalCurrMatchedPix, outlierMask, 0);
	H.release();
}

void PoseEstimation::mapping(string filename, IplImage* image)
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
	int numPoints = curr_x.size();
	int counter = 0;
	CvMat *pointMat = cvCreateMat (3, 1, CV_32FC1);
	CvMat *outPointMat = cvCreateMat (3, 1, CV_32FC1);

	//Opposite globalT
	temp_T->data.fl[0] = -1*currGlobal_T->data.fl[0];
	temp_T->data.fl[1] = -1*currGlobal_T->data.fl[1];
	temp_T->data.fl[2] = -1*currGlobal_T->data.fl[2];

	//R Transpose (opposite globalR)
	temp_R->data.fl[0] = currGlobal_R->data.fl[0];
	temp_R->data.fl[1] = currGlobal_R->data.fl[3];
	temp_R->data.fl[2] = currGlobal_R->data.fl[6];
	temp_R->data.fl[3] = currGlobal_R->data.fl[1];
	temp_R->data.fl[4] = currGlobal_R->data.fl[4];
	temp_R->data.fl[5] = currGlobal_R->data.fl[7];
	temp_R->data.fl[6] = currGlobal_R->data.fl[2];
	temp_R->data.fl[7] = currGlobal_R->data.fl[5];
	temp_R->data.fl[8] = currGlobal_R->data.fl[8];

	for (int i = 0; i < numPoints; i++)
	{
		if(curr_z[i] > 0)
		{
			cvSetReal2D(pointMat, 0, 0, curr_x[i]);
			cvSetReal2D(pointMat, 1, 0, curr_y[i]);
			cvSetReal2D(pointMat, 2, 0, curr_z[i]);
			cvMatMulAdd(temp_R, pointMat, temp_T, outPointMat);
			if(numChannels == 1)
				ss << outPointMat->data.fl[0] << " " << outPointMat->data.fl[1] << " " << outPointMat->data.fl[2] << " " << int(data[i]) << " " << int(data[i]) << " " << int(data[i]) << '\n';
			else if(numChannels == 3)
			{
				ss << outPointMat->data.fl[0] << " " << outPointMat->data.fl[1] << " " << outPointMat->data.fl[2] << " " << int(rData[i]) << " " << int(gData[i]) << " " << int(bData[i]) << '\n';
			}
			counter++;

		}
	}

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

void PoseEstimation::copyGlobal_R_T()
{
	cvCopy(currGlobal_R, prevGlobal_R);
	cvCopy(currGlobal_T, prevGlobal_T);
}

void PoseEstimation::resetPrevPos()
{
	prev_x.clear();
	prev_y.clear();
	prev_z.clear();

	for (int k = 0; k < curr_x.size(); k++)
	{
		prev_x.push_back(curr_x[k]);
		prev_y.push_back(curr_y[k]);
		prev_z.push_back(curr_z[k]);
	}
}

void PoseEstimation::depthToPointCloud(IplImage *depthImage, cv::Mat camIntrinsicMatrix)
{
	//compute point cloud from kinect depth - START
	const int GAMMASIZE =  2048;
	float gamma[GAMMASIZE];
	const float k1 = 1.1863;
	const float k2 = 2842.5;
	const float k3 = 0.1236;

	for (size_t i = 0; i < GAMMASIZE; i++)
		gamma[i] = k3 * tan(i/k2 + k1);


	// camera intrinsic parameters, representative values, see http://nicolas.burrus.name/index.php/Research/KinectCalibration for more info
	float cx = camIntrinsicMatrix.at<double>(0,2);//331.91; //center of projection
	float cy = camIntrinsicMatrix.at<double>(1,2);//259.48; //center of projection
	float fx = camIntrinsicMatrix.at<double>(0,0);//529.24; //focal length in pixels
	float fy = camIntrinsicMatrix.at<double>(1,1);//528.35; //focal length in pixels

	unsigned char *depthData = (unsigned char*)(depthImage->imageData);
	for (int i = 0; i < depthImage->height; i++)
	{
		for (int j = 0; j < depthImage->width; j++)
		{
			float gamma_depth_data = gamma[((short*)depthData)[i*depthImage->width+j]];
			int index = i*depthImage->width+j;
			curr_x[index] = (j - cx) * gamma_depth_data / fx;
			curr_y[index] = (i - cy) * gamma_depth_data / fy;
			curr_z[index] = gamma_depth_data;
		}
	}
}

void PoseEstimation::dispToPointCloud(cv::Mat dispMat,  cv::Mat Q)
{
	cv::Mat xyz;
	reprojectImageTo3D(dispMat, xyz, Q, true);

	int index = 0;
	for (int i = 0; i < dispMat.rows; i++)
	{
		for (int j = 0; j < dispMat.cols; j++)
		{
			cv::Vec3f point = xyz.at<cv::Vec3f>(i, j);
			curr_x[index] = 0.001*point[0];
			curr_y[index] = 0.001*point[1];
			curr_z[index] = -0.001*point[2];
			index++;
		}
	}
	Q.release();
}

void PoseEstimation::readDepthFiles(string filename, vector<string>& depthFiles)
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

void PoseEstimation::savePointCloud(vector<string> filenames, int iteration, IplImage* image)
{
	ifstream fin;
	string filename;
	int index = 0;
	float x, y, z;
	int width = image->width;
	int xPos, yPos;
	stringstream ss;
	ss<<iteration;
	filename = filenames[iteration];
	fin.open(filename.c_str());

	for(int i=0; i<image->width*image->height; i++)
	{
		curr_x[i] = 0;
		curr_y[i] = 0;
		curr_z[i] = 0;
	}

	while(fin.good())
	{
		fin >> yPos;
		fin >> xPos;
		fin >> x;
		fin >> y;
		fin >> z;

		index = yPos*width + xPos;

		curr_x[index] = x;
		curr_y[index] = y;
		curr_z[index] = z;
	}

	fin.close();
}

void PoseEstimation::resizeXYZVectors(int width, int height)
{
	curr_x.resize(width*height);
	curr_y.resize(width*height);
	curr_z.resize(width*height);

	prev_x.resize(width*height);
	prev_y.resize(width*height);
	prev_z.resize(width*height);
}


void PoseEstimation::nearestNeighborMatching()
{
	int numMatches;
	int size;

	cv::flann::Index treeFlannIndex(globalCurrDescriptors, cv::flann::KDTreeIndexParams());

	int k=numNN; // find the 2 nearest neighbors
	cv::Mat results(globalPrevDescriptors.rows, k, CV_32SC1);
	cv::Mat dists(globalPrevDescriptors.rows, k, CV_32FC1);
	treeFlannIndex.knnSearch(globalPrevDescriptors, results, dists, k, cv::flann::SearchParams(flannCheck) ); 

	// Find correspondences by NNDR (Nearest Neighbor Distance Ratio)
	for(int i=0; i < globalPrevKeyPoints.size(); ++i)
	{
		// Apply NNDR
		if(dists.at<float>(i,0) <= nndrRatio * dists.at<float>(i,1))
		{
			//prev frame pixel coordinates for selected matches with prevFrame
			globalPrevMatchedPix.push_back(globalPrevKeyPoints.at(i).pt);
						  
			//curr pixel coordinates for selected matches with prevFrame
			globalCurrMatchedPix.push_back(globalCurrKeyPoints.at(results.at<int>(i,0)).pt);
		}
	}

	results.release();
	dists.release();
}

void PoseEstimation::collect3DMatches(int width)
{
	int currIndex, prevIndex;
	clearMatchedPoints();
	for (int jj=0; jj<globalCurrMatchedPix.size(); jj++)
	{
		//curr frame points
		currIndex = find_curr_index(jj, width);
		
		//previous frame points
		prevIndex = find_prev_index(jj, width);
		
		//fill in the 3D point information
		push_prev3Dpt(prevIndex);
		push_curr3Dpt(currIndex);
	}
}

void PoseEstimation::process(int depthInfo, IplImage* image)
{
	int relativePoseError = 1;
	int i;
	vector<cv::Point3f> currPts, prevPts;
	vector<cv::Point2f> currPix, prevPix;
	float normXDist = 0.0, normYDist = 0.0, normZDist = 0.0;
	float stdXDev = 0.0, stdYDev = 0.0, stdZDev = 0.0;
	bool xTrue, yTrue;
	int numX = 0, numY = 0;
	int counter = 0;

	if(globalCurrDescriptors.type() !=CV_32F)
	{
		globalCurrDescriptors.convertTo(globalCurrDescriptors, CV_32F);
	} 
	if(globalPrevDescriptors.type() != CV_32F)
	{
		globalPrevDescriptors.convertTo(globalPrevDescriptors, CV_32F);
	}

	//Nearest Neighbor Matching
	nearestNeighborMatching();

	//Collect 3D Point Matches
	if(globalCurrMatchedPix.size() >= nbMatches)
		collect3DMatches(image->width);

	matchWeights.clear();

	if(depthInfo != NO_DEPTH)
	{
		//Calculate Norm
		for(i=0; i<globalCurrMatchedPts.size(); i++)
		{
			normXDist += globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x;
			normYDist += globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y;
			normZDist += globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z;
		}
		normXDist/=globalCurrMatchedPts.size();
		normYDist/=globalCurrMatchedPts.size();
		normZDist/=globalCurrMatchedPts.size();

		//Calculate Standard Deviation
		for(i=0; i<globalCurrMatchedPts.size(); i++)
		{
			stdXDev += pow((globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x-normXDist), 2);
			stdYDev += pow((globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y-normYDist), 2);
			stdZDev += pow((globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z-normZDist), 2);
		}
		stdXDev/=globalCurrMatchedPix.size();
		stdYDev/=globalCurrMatchedPix.size();
		stdZDev/=globalCurrMatchedPix.size();

		stdXDev = sqrt(stdXDev);
		stdYDev = sqrt(stdYDev);
		stdZDev = sqrt(stdZDev);

		cout << "Matches Before STD DEV: " << globalCurrMatchedPts.size() << endl;

		//Remove those not within 1 STD DEV of the NORM
		for(i=0; i<globalCurrMatchedPix.size(); i++)
		{
			if((fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x-normXDist) < stdXDev) && (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y-normYDist) < stdYDev) && (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z-normZDist) < stdZDev))
			{
				currPts.push_back(globalCurrMatchedPts[i]);
				currPix.push_back(globalCurrMatchedPix[i]);
				prevPts.push_back(globalPrevMatchedPts[i]);
				prevPix.push_back(globalPrevMatchedPix[i]);

				if(globalCurrMatchedPix[i].y > weightedPortion*(image->height) || globalPrevMatchedPix[i].y > weightedPortion*(image->height))
				{
					matchWeights.push_back(DOUBLE_WEIGHT);
				}
				else
					matchWeights.push_back(SINGLE_WEIGHT);
			}
		}

		globalCurrMatchedPts.clear();
		globalPrevMatchedPts.clear();
		globalCurrMatchedPix.clear();
		globalPrevMatchedPix.clear();
		globalCurrMatchedPts = currPts;
		globalPrevMatchedPts = prevPts;
		globalCurrMatchedPix = currPix;
		globalPrevMatchedPix = prevPix;
		currPts.clear();
		prevPts.clear();
		currPix.clear();
		prevPix.clear();
	}
	else
	{
		//Calculate Norm
		for(i=0; i<globalCurrMatchedPix.size(); i++)
		{
			normXDist += globalCurrMatchedPix[i].x-globalPrevMatchedPix[i].x;
			normYDist += globalCurrMatchedPix[i].y-globalPrevMatchedPix[i].y;
		}
		normXDist/=globalCurrMatchedPix.size();
		normYDist/=globalCurrMatchedPix.size();

		//Calculate Standard Deviation
		for(i=0; i<globalCurrMatchedPts.size(); i++)
		{
			stdXDev += pow((globalCurrMatchedPix[i].x-globalPrevMatchedPix[i].x-normXDist), 2);
			stdYDev += pow((globalCurrMatchedPix[i].y-globalPrevMatchedPix[i].y-normYDist), 2);
		}
		stdXDev/=globalCurrMatchedPix.size();
		stdYDev/=globalCurrMatchedPix.size();

		stdXDev = sqrt(stdXDev);
		stdYDev = sqrt(stdYDev);

		cout << "Matches Before STD DEV: " << globalCurrMatchedPix.size() << endl;

		//Remove those not within 1 STD DEV of the NORM
		for(i=0; i<globalCurrMatchedPix.size(); i++)
		{
			if((fabs(globalCurrMatchedPix[i].x-globalPrevMatchedPix[i].x-normXDist) < stdXDev) && (fabs(globalCurrMatchedPix[i].y-globalPrevMatchedPix[i].y-normYDist) < stdYDev))
			{
				currPts.push_back(globalCurrMatchedPts[i]);
				currPix.push_back(globalCurrMatchedPix[i]);
				prevPts.push_back(globalPrevMatchedPts[i]);
				prevPix.push_back(globalPrevMatchedPix[i]);

				if(globalCurrMatchedPix[i].y > weightedPortion*(image->height) || globalPrevMatchedPix[i].y > weightedPortion*(image->height))
				{
					matchWeights.push_back(DOUBLE_WEIGHT);
				}
				else
					matchWeights.push_back(SINGLE_WEIGHT);
			}
		}

		globalCurrMatchedPts.clear();
		globalPrevMatchedPts.clear();
		globalCurrMatchedPix.clear();
		globalPrevMatchedPix.clear();
		globalCurrMatchedPts = currPts;
		globalPrevMatchedPts = prevPts;
		globalCurrMatchedPix = currPix;
		globalPrevMatchedPix = prevPix;
		currPts.clear();
		prevPts.clear();
		currPix.clear();
		prevPix.clear();
	}

	cout << "Matches After STD DEV: " << globalCurrMatchedPts.size() << endl;

	if (globalCurrMatchedPix.size() >= nbMatches)
	{
		find_homography();

		if(depthInfo == NO_DEPTH)
		{
			cout << endl;
			cout << "*****************************************************************" << endl;
			cout << "* WARNING: Trying to use depth info for estimateRelativePose... *" << endl;
			cout << "* Continuing on without pose information...                     *" << endl;
			cout << "*****************************************************************" << endl;
			cout << endl;

			relativePoseError = 1;
		}
		else
		{
			relativePoseError = estimateRelativePose();
		}
	}
	else
	{
		cout << "~~~~~~~~~~ NO POSE ESTIMATION! NOT ENOUGH GOOD POINTS ~~~~~~~~~~" << endl;
	}

	if (relativePoseError == 0)
	{
		composePose();
		printCurrentGlobal_R_T();
		mapping(string("results/PointCloud.txt"), image);
		copyGlobal_R_T();
	}
}

void PoseEstimation::clear()
{
	currMatchedPoints.clear();
	prevMatchedPoints.clear();
	prevMatchedPixels.clear();
	currMatchedPixels.clear();
	prevObjectData.release();
	currObjectData.release();
	prevKeypoints.clear();
	outlierMask.release();
}


int PoseEstimation::searchPointProj(Point2f find)
{
	int match;
	int end;
	for(int i=0; i<projectedPoints.size(); i++)
	{
		match = i;
		end = projectedPoints[i].pixels.size()-1;

		if(fabs(projectedPoints[i].pixels[end].x-find.x)+fabs(projectedPoints[i].pixels[end].y-find.y)<0.001)
			return match;
	}
	return -1;
}

void PoseEstimation::appendPointProj(Point2f pix, Point3f pt, int loc, int frameIndex)
{
	projectedPoints[loc].pixels.push_back(pix);
	projectedPoints[loc].allPoints.push_back(pt);
	projectedPoints[loc].frames.push_back(frameIndex);
}

void PoseEstimation::savePointProj(int frameIndex, int firstFrame)
{
	int index, match;
	pointProjection tempPoints;

	for(int i=0; i<validPix.size(); i++)
	{
		index = validPix[i];
		tempPoints.point = globalCurrMatchedPts[index];
		tempPoints.allPoints.push_back(globalPrevMatchedPts[index]);
		tempPoints.pixels.push_back(globalPrevMatchedPix[index]);
		tempPoints.pixels.push_back(globalCurrMatchedPix[index]);
		tempPoints.frames.push_back(frameIndex);
		tempPoints.frames.push_back(frameIndex+1);

		if(frameIndex == firstFrame)
			projectedPoints.push_back(tempPoints);
		else
		{
			match = searchPointProj(globalPrevMatchedPix[index]);
			if(match != -1)
			{
				appendPointProj(globalCurrMatchedPix[index], globalCurrMatchedPts[index], match, frameIndex+1);
			}
			else
			{
				projectedPoints.push_back(tempPoints);
			}
		}
		tempPoints.allPoints.clear();
		tempPoints.pixels.clear();
		tempPoints.frames.clear();
	}
}

void PoseEstimation::writePointProj(string pointProjFile)
{
	ofstream fout;

	fout.open(pointProjFile.c_str());
	for(int i=0; i<projectedPoints.size(); i++)
	{
		fout << projectedPoints[i].point.x << " " << projectedPoints[i].point.y << " "<< projectedPoints[i].point.z << " " << projectedPoints[i].frames.size() << " ";
		for(int j=0; j<projectedPoints[i].frames.size(); j++)
		{
			fout << projectedPoints[i].frames[j] << " " << projectedPoints[i].pixels[j].x << " " << projectedPoints[i].pixels[j].y << " ";
		}
		fout << endl;
	}
	fout.close();
}


