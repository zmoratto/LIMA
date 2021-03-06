/*weights.h: header for weights program */
#ifndef WEIGHTS_H
#define WEIGHTS_H
//the header goes here


//includes
#include "lidar_tracks/tracks.h"
#include "featuresLOLA.h"
#include <vector>
#include <stdio.h>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

//function prototypes
int MakeLinearWeights( vector< LOLAShot > & trackPts, const int &halfWindow);
void SaveWeights(vector< LOLAShot>& trackPts, string filename);
void ComputeWeights( vector< vector<LOLAShot> > & track_pts, int halfWindow, float topPercent,  string filename);
void ResetWeights( vector< vector<LOLAShot> >& trackPts );

#endif 
