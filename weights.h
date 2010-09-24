/*weights.h: header for weights program */
#ifndef WEIGHTS_H_GUARD
#define WEIGHTS_H_GUARD
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
int simple_function_test(int N);
int correct_invalid_linear(vector<vector< LOLAShot> > & track_pts);
int check_acceptable_computes(vector< vector<LOLAShot > > & track_pts,int edge, int & num_valid);
int grad_filter(vector< vector<LOLAShot > > & track_pts, int edge);
int make_weights_not_war( vector<vector< LOLAShot > > & track_pts, int edge, float take_p, int& num_valid, float& take_thresh);
int write_refl_and_weights( vector< vector< LOLAShot> > & track_pts, string & f_name);



#endif 
