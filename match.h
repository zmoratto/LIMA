// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef MATCH_H
#define MATCH_H

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include "coregister.h"
#include "tracks.h"
#include "weights.h"
#include "util.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

#define DEFAULT_SEARCH_TRANS_WINDOW 20
#define DEFAULT_SEARCH_TRANS_STEP 5.0
#define DEFAULT_SEARCH_THETA_WINDOW (M_PI / 10)
#define DEFAULT_SEARCH_THETA_STEP (M_PI / 40)

Matrix3x3 find_track_homography(vector<vector<AlignedLOLAShot> > aligned, vector<vector<AlignedLOLAShot> > aligned2);
Matrix3x3 find_image_homography(char* image1, char* image2, float lonstart, float lonend, float latstart, float latend);

vector<vector<AlignedLOLAShot> > align_to_image(vector<vector<LOLAShot> > & trackPts, ImageView<PixelGray<float> > cubImage,
		Matrix3x3 & trans, string image_file = "");
vector<vector<AlignedLOLAShot> > align_to_image_pyramid(vector<vector<LOLAShot> > & trackPts, const string & image_file,
		Matrix3x3 & trans, string outputImage = "");

float compute_transform_error(vector<vector<AlignedLOLAShot> > & tracks, int* numpoints=NULL);
void find_track_transforms(vector<vector<AlignedLOLAShot> > & tracks, string cubFile);
Matrix3x3 find_tracks_transform(vector<vector<AlignedLOLAShot> > & tracks, ImageView<PixelGray<float> > & cub,
		int transSearchWindow=DEFAULT_SEARCH_TRANS_WINDOW, int transSearchStep=DEFAULT_SEARCH_TRANS_STEP, 
		float thetaSearchWindow=DEFAULT_SEARCH_THETA_WINDOW, float thetaSearchStep=DEFAULT_SEARCH_THETA_STEP);

Matrix3x3 gauss_newton_homography(vector<AlignedLOLAShot> & track, ImageView<PixelGray<float> > cubImage,
		Matrix3x3 matrix=Matrix3x3(1,0,0,0,1,0,0,0,1));
Matrix3x3 gauss_newton_affine(vector<AlignedLOLAShot> & track, ImageView<PixelGray<float> > cubImage,
		Matrix3x3 matrix=Matrix3x3(1,0,0,0,1,0,0,0,1));

float ComputeScaleFactor(vector<float> allImgPts, vector<float> reflectance);
float ComputeScaleFactor(vector<Vector3> allImgPts, vector<float> reflectance);

void GenerateInitTransforms( vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings);
void GetBestTransform(vector<Vector<float, 6> > &finalTransfArray,  vector<float> &finalMatchingErrorArray, 
                      Vector<float, 6> &optimalTransf, float &optimalError);

//void InitMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename,  
//			       CoregistrationParams coregistrationParams,  
//			       vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransf, 
//			       vector<float> &matchingError);

vector<Vector4> FindMatches2D(vector<vector<AlignedLOLAShot> > &trackPts, string cubFilename, 
                              Vector2 matchWindowHalfSize, int numSamples, vector<float> &errorArray);
void EstimateAffineTransform(vector<Vector<float, 6> > &initTransfArray, vector<vector<LOLAShot> > &trackPts, 
                             vector<Vector4> matchArray, vector<float> errorArray);

//void UpdateMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename,  
//			         int numMaxIter, vector<Vector<float, 6> >initTransfArray, vector<float> &initErrorArray, 
//                               vector<Vector<float, 6> >&finalTransfArray, vector<float> &errorArray, Vector2 &centroid );

void EstimateMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename, Vector2 matchWindowHalfSize, 
                   vector<Vector<float, 6> > &initTransfArray,  vector<float> &matchingErrorArray);

float ComputeMatchingError(vector<float> reflectancePts, vector<float>imgPts);


void SaveReportFile(vector<vector<LOLAShot> > &trackPts, vector<Vector<float, 6> >finalTransfArray,  
                    vector<float> errorArray, string matchResultsFilename);

/*
void InitMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
			ModelParams modelParams, CoregistrationParams coregistrationParams,  
			vector<Vector<float, 6> >initTransfArray, Vector<float, 6> &finalTransf, 
			float *matchingError);
*/
/*
void UpdateMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename, 
                     ModelParams modelParams,  int numMaxIter, 
                     vector<Vector<float, 6> > initTransfArray,  vector<Vector<float, 6> >&finalTransfArray, 
                     vector<float> &errorArray);
*/
/*
void UpdateMatchingParamsLIMA_MP(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
				 ModelParams modelParams, CoregistrationParams coregistrationParams,
				 vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
				 vector<float> &errorArray );
*/
/*
void UpdateMatchingParamsLIDEM_MP(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
			          ModelParams modelParams, CoregistrationParams coregistrationParams,
			          vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
				  vector<float> &errorArray );
*/
#endif//MATCH_H
