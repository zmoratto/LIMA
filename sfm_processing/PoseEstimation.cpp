#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "opencv2/calib3d/calib3d.hpp"
#include <pcl/io/pcd_io.h>
#include "pcl/common/common_headers.h"
#include <pcl/correspondence.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/registration/transformation_estimation_lm.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>

#include <Eigen/Dense>

#include "PoseEstimation.h"

using namespace Eigen;
using namespace std;

PoseEstimation::PoseEstimation()
{
	//Clear Point Vectors
	currMatchedPoints.clear();
	prevMatchedPoints.clear();
	prevMatchedPixels.clear();
	currMatchedPixels.clear();
	prevKeypoints.clear();
	globalCurrMatchedPts.clear();
	globalPrevMatchedPts.clear();
	globalPrevMatchedPix.clear();
	globalCurrMatchedPix.clear();
	globalCurrKeyPoints.clear();
	globalPrevKeyPoints.clear();

	//Previous Global R and T With First Frame
	prevGlobal_R = cvCreateMat(3,3,CV_32F);
	prevGlobal_T = cvCreateMat(3,1,CV_32F);
	//R
	cvSetReal2D(prevGlobal_R, 0, 0, 1.0);
	cvSetReal2D(prevGlobal_R, 0, 1, 0.0);
	cvSetReal2D(prevGlobal_R, 0, 2, 0.0);
	cvSetReal2D(prevGlobal_R, 1, 0, 0.0);
	cvSetReal2D(prevGlobal_R, 1, 1, 1.0);
	cvSetReal2D(prevGlobal_R, 1, 2, 0.0);
	cvSetReal2D(prevGlobal_R, 2, 0, 0.0);
	cvSetReal2D(prevGlobal_R, 2, 1, 0.0);
	cvSetReal2D(prevGlobal_R, 2, 2, 1.0);
	//T
	cvSetReal2D(prevGlobal_T, 0, 0, 0.0);
	cvSetReal2D(prevGlobal_T, 1, 0, 0.0);
	cvSetReal2D(prevGlobal_T, 2, 0, 0.0);	

	//Current Global R and T
	currGlobal_R = cvCreateMat(3,3,CV_32F);
	currGlobal_T = cvCreateMat(3,1,CV_32F);

	//Relative R and T
	relative_R = cvCreateMat (3, 3, CV_32FC1);
	relative_T = cvCreateMat(3,1,CV_32F);

	//Clear All Other Vectors
	projectedPoints.clear();
	validPix.clear();
	matchWeights.clear();
	curr_x.clear();
	curr_y.clear();
	curr_z.clear();
	prev_x.clear();
	prev_y.clear();
	prev_z.clear();
}

PoseEstimation::~PoseEstimation()
{
	//Release all Mats
	outlierMask.release();
	cvReleaseMat(&prevGlobal_R);
	cvReleaseMat(&prevGlobal_T);
	cvReleaseMat(&currGlobal_R);
	cvReleaseMat(&currGlobal_T);
	cvReleaseMat(&relative_R);
	cvReleaseMat(&relative_T);

	//Release all Vectors
	currMatchedPoints.clear();
	prevMatchedPoints.clear();
	prevMatchedPixels.clear();
	currMatchedPixels.clear();
	prevKeypoints.clear();
	globalCurrMatchedPts.clear();
	globalPrevMatchedPts.clear();
	globalPrevMatchedPix.clear();
	globalCurrMatchedPix.clear();
	globalCurrKeyPoints.clear();
	globalPrevKeyPoints.clear();
	projectedPoints.clear();
	validPix.clear();
	matchWeights.clear();
	curr_x.clear();
	curr_y.clear();
	curr_z.clear();
	prev_x.clear();
	prev_y.clear();
	prev_z.clear();
}

//Estimates relative R and T between two point cloud sets. These point cloud matches are in globalMatchedPts.
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
	CvMat* inMat = cvCreateMat(3,1,CV_32FC1);
	CvMat* outMat = cvCreateMat(3,1, CV_32FC1);

	//Find Translation, Previous Center, and Current Center of Point Clouds
	for (i = 0; i < numMatches; i++)
	{
		//If not an outlier and the matching coordinates are within given thresholds of each other, they are a valid match.
		if ((outlierMask.at<bool>(i) == 1) && (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z) < zDist) && 
                    (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y) < yDist) && (fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x) < xDist))
		{
			//Weight the translation matrix by where the match is located. More weight if the match is closer to the camera.
			if(matchWeights[i] != NO_WEIGHT)
			{
				currCenter = currCenter + (globalCurrMatchedPts[i]);
				prevCenter = prevCenter + (globalPrevMatchedPts[i]);
				validPix.push_back(i);
				numValidMatches++;
			}
			numWeightedMatches+=matchWeights[i]; //Collect the weighted values of the matches
		}
	}

	cout << "NumMatches: " << numMatches << endl; //Total number of matches
	cout << "NumValidMatches: " << numValidMatches << endl; //Number of matches that were used to calculate T and the centers.
	cout << "NumWeights: " << numWeightedMatches << endl; //Collected weights of the matches used.

	//Normalize the translation and centers.
	if (numValidMatches > 0)
	{
		currCenter.x = currCenter.x/numValidMatches;
		currCenter.y = currCenter.y/numValidMatches;
		currCenter.z = currCenter.z/numValidMatches;
		prevCenter.x = prevCenter.x/numValidMatches;
		prevCenter.y = prevCenter.y/numValidMatches;
		prevCenter.z = prevCenter.z/numValidMatches;

	}
	else
	{
		return 1; //Return ERROR, not enough matches.
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

	//Compute the rotation matrix Relative_R
	for (i = 0; i < numMatches; i++)
	{
		//If the match is a valid match.
		if ((outlierMask.at<bool>(i) == 1) && (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z) < zDist) && 
                    (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y) < yDist) && (fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x) < xDist))
		{
			//Subtract the centers from each matched point and compute the A matrix.
			if(matchWeights[i] != NO_WEIGHT)
			{
				a00 = a00 + (globalCurrMatchedPts[i].x - currCenter.x)*(globalPrevMatchedPts[i].x - prevCenter.x);
				a10 = a10 + (globalCurrMatchedPts[i].y - currCenter.y)*(globalPrevMatchedPts[i].x - prevCenter.x);
				a20 = a20 + (globalCurrMatchedPts[i].z - currCenter.z)*(globalPrevMatchedPts[i].x - prevCenter.x);

				a01 = a01 + (globalCurrMatchedPts[i].x - currCenter.x)*(globalPrevMatchedPts[i].y - prevCenter.y);
				a11 = a11 + (globalCurrMatchedPts[i].y - currCenter.y)*(globalPrevMatchedPts[i].y - prevCenter.y);
				a21 = a21 + (globalCurrMatchedPts[i].z - currCenter.z)*(globalPrevMatchedPts[i].y - prevCenter.y);

				a02 = a02 + (globalCurrMatchedPts[i].x - currCenter.x)*(globalPrevMatchedPts[i].z - prevCenter.z);
				a12 = a12 + (globalCurrMatchedPts[i].y - currCenter.y)*(globalPrevMatchedPts[i].z - prevCenter.z);
				a22 = a22 + (globalCurrMatchedPts[i].z - currCenter.z)*(globalPrevMatchedPts[i].z - prevCenter.z);
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

	//If we have enough valid matches, we can compute SVD of A.
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
			return 1; //ERROR, determinant == 0.
		}
		else //Use Kabsch Algorithm for determining the rotation matrix between two sets of point clouds.
		{
			if (det > 0)
			{
				//Use the decomposition of A to compute the relative_R matrix
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

		//Find Translation New Way
		for (i = 0; i < numMatches; i++)
		{
			if ((outlierMask.at<bool>(i) == 1) && (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z) < zDist) && (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y) < yDist) && (fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x) < xDist))
			{
				cvSetReal2D(inMat, 0, 0, globalPrevMatchedPts[i].x);
				cvSetReal2D(inMat, 1, 0, globalPrevMatchedPts[i].y);
				cvSetReal2D(inMat, 2, 0, globalPrevMatchedPts[i].z);
				cvMatMul(relative_R, inMat, outMat);
				translation.x += matchWeights[i]*(globalCurrMatchedPts[i].x - outMat->data.fl[0]);
				translation.y += matchWeights[i]*(globalCurrMatchedPts[i].y - outMat->data.fl[1]);
				translation.z += matchWeights[i]*(globalCurrMatchedPts[i].z - outMat->data.fl[2]);
			}
		}

		translation.x = translation.x/numWeightedMatches;
		translation.y = translation.y/numWeightedMatches;
		translation.z = translation.z/numWeightedMatches;
		cvSetReal2D(relative_T, 0, 0, translation.x);
		cvSetReal2D(relative_T, 1, 0, translation.y);
		cvSetReal2D(relative_T, 2, 0, translation.z);

		return 0;
	}
	else
	{
		cvReleaseMat(&AMat);
		return 1; //ERROR, not enough valid matches
	}
}

//Uses relative R and T and prevGlobal R and T to compute the currGlobal R and T with respect to the first frame.
void PoseEstimation::composePose()
{
	//Compute Global R and Global T using relative R and relative T.
	cvMatMul(prevGlobal_R, relative_R, currGlobal_R);
	cvMatMulAdd(prevGlobal_R, relative_T, prevGlobal_T, currGlobal_T);
}

//Wrapper for OpenCV findHomography and findFundamentalMat for two sets of pixels
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

//Copy global R and T to previous Global R and T
void PoseEstimation::copyGlobal_R_T()
{
	cvCopy(currGlobal_R, prevGlobal_R);
	cvCopy(currGlobal_T, prevGlobal_T);
}

//Copy the current 3D points to the previous 3D points and clean out current points.
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
		curr_x[k] = 0.0;
		curr_y[k] = 0.0;
		curr_z[k] = 0.0;
	}
}

//Allocate memory for 3D point holders
void PoseEstimation::allocatePCMemory(int width, int height)
{
	curr_x.resize(width*height);
	curr_y.resize(width*height);
	curr_z.resize(width*height);

	prev_x.resize(width*height);
	prev_y.resize(width*height);
	prev_z.resize(width*height);
}

//Flann Nearest Neighbor Matching
void PoseEstimation::nearestNeighborMatching(IplImage* image)
{
	int numMatches;
	int size;
	cv::KeyPoint pPt, cPt;
	float cX, cY, pX, pY;
	bool match = false;
	int k=numNN; // find the 2 nearest neighbors

	//NNDR Matching
	if(matchingMethod == 0)
	{
		cv::flann::Index treeFlannIndex(globalCurrDescriptors, cv::flann::KDTreeIndexParams());
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
	else //Marrying Matches
	{
		//Backward
		cv::flann::Index backwardIndex(globalCurrDescriptors, cv::flann::KDTreeIndexParams());

		cv::Mat resultsBack(globalPrevDescriptors.rows, 1, CV_32SC1);
		cv::Mat distsBack(globalPrevDescriptors.rows, 1, CV_32FC1);
		backwardIndex.knnSearch(globalPrevDescriptors, resultsBack, distsBack, 1, cv::flann::SearchParams(flannCheck)); 

		//Forward
		cv::flann::Index forwardIndex(globalPrevDescriptors, cv::flann::KDTreeIndexParams());
		cv::Mat resultsForward(globalCurrDescriptors.rows, 1, CV_32SC1);
		cv::Mat distsForward(globalCurrDescriptors.rows, 1, CV_32FC1);
		forwardIndex.knnSearch(globalCurrDescriptors, resultsForward, distsForward, 1, cv::flann::SearchParams(flannCheck));
	
		// Find correspondences by NNDR (Nearest Neighbor Distance Ratio)
		for(int i=0; i < globalPrevKeyPoints.size(); ++i)
		{
			pPt = globalPrevKeyPoints.at(i);
			pX = pPt.pt.x;
			pY = pPt.pt.y;

			cPt = globalCurrKeyPoints.at(resultsBack.at<int>(i,0));
			cX = cPt.pt.x;
			cY = cPt.pt.y;

			for(int j=0; j<globalCurrKeyPoints.size(); j++)
			{
				if((globalCurrKeyPoints.at(j).pt.x == cX) && (globalCurrKeyPoints.at(j).pt.y == cY) && (globalPrevKeyPoints.at(resultsForward.at<int>(j,0)).pt.x == pX) && (globalPrevKeyPoints.at(resultsForward.at<int>(j,0)).pt.y == pY))
				{
					globalPrevMatchedPix.push_back(globalPrevKeyPoints.at(i).pt);
					globalCurrMatchedPix.push_back(globalCurrKeyPoints.at(resultsBack.at<int>(i,0)).pt);
					break;
				}
			}
		}

		resultsBack.release();
		distsBack.release();
		resultsForward.release();
		distsForward.release();
	}

	//Collect 3D Point Matches
	if(globalCurrMatchedPix.size() >= nbMatches)
		collect3DMatches(image->width);
}

//Using Pixel Matches, find corresponding 3D Matches
void PoseEstimation::collect3DMatches(int width)
{
	int currIndex, prevIndex;
	currMatchedPoints.clear();
	prevMatchedPoints.clear();

	for (int jj=0; jj<globalCurrMatchedPix.size(); jj++)
	{
		//curr frame points
		currIndex = int(globalCurrMatchedPix[jj].x)+int(globalCurrMatchedPix[jj].y)*width;//find_curr_index(jj, width);
		
		//previous frame points
		prevIndex = int(globalPrevMatchedPix[jj].x)+int(globalPrevMatchedPix[jj].y)*width;//find_prev_index(jj, width);
		
		//fill in the 3D point information
		cv::Point3f pt3;
		pt3.x = prev_x[prevIndex];
		pt3.y = prev_y[prevIndex];  
		pt3.z = prev_z[prevIndex];
		globalPrevMatchedPts.push_back(pt3);

		pt3.x = curr_x[currIndex];
		pt3.y = curr_y[currIndex];  
		pt3.z = curr_z[currIndex];
		globalCurrMatchedPts.push_back(pt3);
	}
}

//Perform Point Cloud Alignment on matches from two images
void PoseEstimation::process(int depthInfo, IplImage* image)
{
	int relativePoseError = 1;
	int i;
	cv::Mat currTemp, prevTemp;
	clock_t t1, t2;

	//Nearest Neighbor Matching
	nearestNeighborMatching(image);

	//remove outliers
	removeOutliers(depthInfo, image);

	if (globalCurrMatchedPix.size() >= nbMatches)
	{
		find_homography();

		if(depthInfo == NO_DEPTH)
			relativePoseError = 1;
		else
			relativePoseError = estimateRelativePose();
	}
	else
		cout << "* NO POSE ESTIMATION! NOT ENOUGH GOOD POINTS *" << endl;

	if (relativePoseError == 0)
	{
		composePose();
		copyGlobal_R_T();
	}
}

//Clear variables for processing of new tile
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

//Search for point projection matches
int PoseEstimation::searchPointProj(cv::Point2f find)
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

//Add new point projection
void PoseEstimation::appendPointProj(cv::Point2f pix, cv::Point3f pt, int loc, int frameIndex)
{
	projectedPoints[loc].pixels.push_back(pix);
	projectedPoints[loc].allPoints.push_back(pt);
	projectedPoints[loc].frames.push_back(frameIndex);
}

//Save All new points into point projection
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
		tempPoints.frames.push_back(frameIndex-1);
		tempPoints.frames.push_back(frameIndex);

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

//Write projected points to a file
void PoseEstimation::writePointProj(string& pointProjFile)
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

//Remove duplicate key points before matching
void PoseEstimation::removeDuplicates(IplImage* image)
{
	CvMat* im = cvCreateMat(image->height,image->width,CV_32F);
	vector<int> repeats;
	vector<cv::KeyPoint> tempKeyPoints;
	int x, y, index, width = image->width;
	tempKeyPoints.clear();
	repeats.clear();
	cv::Mat tempMat;

	for(int i=0; i<(image->width)*(image->height); i++)
		im->data.fl[i] = 0.0;

	for(int i=0; i<globalCurrKeyPoints.size(); i++)
	{
		x = globalCurrKeyPoints[i].pt.x;
		y = globalCurrKeyPoints[i].pt.y;

		index = y*width + x;

		if(im->data.fl[index] != 255)
		{
			tempKeyPoints.push_back(globalCurrKeyPoints[i]);
			im->data.fl[index] = 255;
		}
		else
			repeats.push_back(i);
	}

	tempMat.create(tempKeyPoints.size(), globalCurrDescriptors.cols, CV_32F);

	int j=0, cnt=0;
	for(int i=0; i<globalCurrKeyPoints.size(); i++)
	{
		if(repeats[j] != i)
		{
			cv::Mat dest(tempMat.rowRange(cnt, cnt+1));
			globalCurrDescriptors.row(i).copyTo(dest);
			cnt++;
			dest.release();
		}
		else
			j++;
	}

	globalCurrDescriptors.release();
	globalCurrDescriptors = tempMat.clone();

	globalCurrDescriptors.convertTo(globalCurrDescriptors, CV_32F);

	globalCurrKeyPoints.clear();
	globalCurrKeyPoints = tempKeyPoints;

	tempKeyPoints.clear();
	tempMat.release();

	repeats.clear();
	cvReleaseMat(&im);
}

//Remove outliers after initial nearest neighbor matching
void PoseEstimation::removeOutliers(int depthInfo, IplImage* image)
{
	vector<cv::Point3f> currPts, prevPts;
	vector<cv::Point2f> currPix, prevPix;
	float normXDist = 0.0, normYDist = 0.0, normZDist = 0.0;
	float stdXDev = 0.0, stdYDev = 0.0, stdZDev = 0.0;
	bool xTrue, yTrue;
	int numX = 0, numY = 0;
	int counter = 0;
	int i;

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
			if((fabs(globalCurrMatchedPts[i].x-globalPrevMatchedPts[i].x-normXDist) < stdXDev) && 
                           (fabs(globalCurrMatchedPts[i].y-globalPrevMatchedPts[i].y-normYDist) < stdYDev) && 
                           (fabs(globalCurrMatchedPts[i].z-globalPrevMatchedPts[i].z-normZDist) < stdZDev))
			{
				currPts.push_back(globalCurrMatchedPts[i]);
				currPix.push_back(globalCurrMatchedPix[i]);
				prevPts.push_back(globalPrevMatchedPts[i]);
				prevPix.push_back(globalPrevMatchedPix[i]);

				//WEIGHTING
				if(globalCurrMatchedPix[i].y > weightedPortion*(image->height) && globalPrevMatchedPix[i].y > weightedPortion*(image->height))
				{
					matchWeights.push_back(DOUBLE_WEIGHT);
				}
				else
					matchWeights.push_back(SINGLE_WEIGHT);
			}
		}
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
			if ((fabs(globalCurrMatchedPix[i].x-globalPrevMatchedPix[i].x-normXDist) < stdXDev) && (fabs(globalCurrMatchedPix[i].y-globalPrevMatchedPix[i].y-normYDist) < stdYDev))
			{
				currPts.push_back(globalCurrMatchedPts[i]);
				currPix.push_back(globalCurrMatchedPix[i]);
				prevPts.push_back(globalPrevMatchedPts[i]);
				prevPix.push_back(globalPrevMatchedPix[i]);
			}
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

