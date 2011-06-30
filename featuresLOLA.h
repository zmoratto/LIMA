/*featuresLOLA.h: header for LOLA features program */
#ifndef FEATURES_LOLA_H
#define FEATURES_LOLA_H
//the header goes here

//includes
#include "tracks.h"
#include <vector>
#include <stdio.h>
#include <string>
#include <iostream>
#include <algorithm>
using namespace std;

//this function loops through each set of track points separatly, 
//if we see the following pattern of points (valid, in-valid, valid)
//the non-valid point will have it's reflectance set to the average of the two valid points

//above non-valid was defined as reflectance = -1 
int InterpolateInvalidPoint(vector< LOLAShot>& trackPts);

//function loops through points setting the LOLAShot integer field .calc_acp to "1" 
//if the point is valid & there are sufficent valid points on either side of that point compute a gradient.  
//otherwise .calc_acp is set to 0 if either the point is invalid or there are insufficient 
//valid points bordering that point to compute a gradient
//halfWindow = (length of the derivative filter-1)/2
//numValid = the number of points for which the derivative filter is computed
 int FindValidPoints(vector<LOLAShot > & trackPts,int halfWindow, int &numValid);

//calculates the gradient of reflectance using a filter of shape: [ -1x--edge-- 0  +1x--edge--]
int ComputeGradient(vector<LOLAShot >& trackPts, int halfWindow);

int ComputeSalientReflectanceFeatures( vector< LOLAShot > & trackPts, float topPercent, int numValid);
void ComputeSalientLOLAFeature(vector<LOLAShot > & trackPts, vector<float> filter, float salientFeatureThresh);


//this function will be grandfathered
void ComputeSalientLOLAFeature(vector<LOLAShot > & trackPts,int halfWindow, float topPercent);
#endif
