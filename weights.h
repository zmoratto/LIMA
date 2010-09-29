/*weights.h: header for weights program */
#ifndef WEIGHTS_H
#define WEIGHTS_H
//the header goes here


//includes
#include "tracks.h"
#include <vector>
#include <stdio.h>
#include <string>
#include <iostream>
#include <algorithm>
using namespace std;

//function prototypes

int InterpolateInvalidPoint(vector<vector< LOLAShot> > & track_pts);
int FindValidPoints(vector< vector<LOLAShot > > & track_pts,int edge, int & num_valid);
int ComputeGradient(vector< vector<LOLAShot > > & track_pts, int edge);
int ComputeSalientFeatures( vector<vector< LOLAShot > > & track_pts, int edge, float take_p, int& num_valid, float& take_thresh);
int SaveWeights( vector< vector< LOLAShot> > & track_pts, string & f_name);

#endif 
