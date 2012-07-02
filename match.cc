// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <omp.h>

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
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>
#include <vw/Math/LinearAlgebra.h>
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

float compute_transform_error(vector<AlignedLOLAShot> & track, int* numpoints=NULL)
{
	float err = 0.0;
	int num_points = 0;
	for (unsigned int j = 0; j < track.size(); j++)
	{
		if (track[j].synth_image == -1 || track[j].image == -1)
			continue;
		err += pow(track[j].synth_image - track[j].image, 2);
		num_points++;
	}
	if (numpoints != NULL)
		*numpoints = num_points;
	if (num_points == 0)
		return -1;

	return err / num_points;
}

float compute_transform_error(vector<vector<AlignedLOLAShot> > & tracks, int* numpoints=NULL)
{
	float err = 0.0;
	int num_points = 0;
	for (unsigned int i = 0; i < tracks.size(); i++)
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			if (tracks[i][j].synth_image == -1 || tracks[i][j].image== -1)
				continue;
			err += pow(tracks[i][j].synth_image - tracks[i][j].image, 2);
			num_points++;
		}
	if (numpoints != NULL)
		*numpoints = num_points;
	if (num_points == 0)
		return -1;

	return err / num_points;
}

Matrix3x3 find_tracks_transform(vector<vector<AlignedLOLAShot> > & tracks, ImageView<PixelGray<float> > & cub,
		int transSearchWindow, int transSearchStep, float thetaSearchWindow, float thetaSearchStep)
{
	float mid_x = cub.cols() / 2;
	float mid_y = cub.rows() / 2;
	Matrix3x3 center(1, 0, -mid_x, 0, 1, -mid_y, 0, 0, 1);
	Matrix3x3 uncenter(1, 0, mid_x, 0, 1, mid_y, 0, 0, 1);

	Matrix3x3 best;
	float best_score = INFINITY;


	for (float tt = -thetaSearchWindow; tt <= thetaSearchWindow; tt += thetaSearchStep)
	{
		Matrix3x3 purerot(cos(tt), -sin(tt), 0, sin(tt), cos(tt), 0, 0, 0, 1);
		Matrix3x3 rot = uncenter * purerot * center;
		for (int xt = -transSearchWindow; xt <= transSearchWindow; xt += transSearchStep)
			for (int yt = -transSearchWindow; yt <= transSearchWindow; yt += transSearchStep)
			{
				Matrix3x3 trans(1, 0, xt, 0, 1, yt, 0, 0, 1);
				transform_tracks(tracks, trans * rot, cub);
				float score = compute_transform_error(tracks);
				if (score < best_score)
				{
					//printf("%g %d %d %g\n", tt, xt, yt, score);
					best_score = score;
					best = trans * rot;
				}
			}
		//printf("%g\n", tt);
	}
	//printf("Best x: %g y: %g theta: %g\n", best_trans_x, best_trans_y, best_rot);
	
	transform_tracks(tracks, best, cub);

	return best;
}

// J_f(B)_{ij} = dF_i / dB_j
// transform may be marked as invalid if points have left the image window
Matrix<double> compute_jacobian(vector<AlignedLOLAShot> & track, Matrix3x3 B, ImageView<PixelGray<float> > img, int num_points)
{
	Matrix<double> J(num_points, 6);
	int index = 0;
	for (unsigned int i = 0; i < track.size(); i++)
	{
		// TODO: we may ignore some points that weren't in the image when num_points was computed,
		// but are with the new transformation
		if (track[i].image == -1 || track[i].synth_image == -1 || index >= num_points)
			continue;

		int x = track[i].image_x;
		int y = track[i].image_y;

		// finite difference with two points
		//double dx = -0.5 * (img(x + 1, y) - img(x - 1, y));
		//double dy = -0.5 * (img(x, y + 1) - img(x, y - 1));
		double h = 1.0;
		
		// fininte difference with three points to right and left
		double dx, dy;
		if (x >= img.rows() - 3 || x <= 3 || y >= img.cols() - 3 || y <= 3)
			dx = 0.0, dy = 0.0; // TODO: something better
		else
		{
			dx = (-1.0/60) * img(x-3, y) + (3.0/20) * img(x-2, y) - (3.0 / 4) * img(x-1, y) +
				(3.0/4) * img(x+1, y) - (3.0 / 20) * img(x+2, y) + (1.0 / 60) * img(x+3, y);
			dy = (-1.0/60) * img(x, y-3) + (3.0/20) * img(x, y-2) - (3.0 / 4) * img(x, y-1) +
				(3.0/4) * img(x, y+1) - (3.0 / 20) * img(x, y+2) + (1.0 / 60) * img(x, y+3);
			dx = -dx;
			dy = -dy;
		}
		
		// compute the derivatives
		J(index, 0) = dx / (h / track[i].imgPt[2].x);
		J(index, 1) = dx / (h / track[i].imgPt[2].y);
		J(index, 2) = dx / h;
		J(index, 3) = dy / (h / track[i].imgPt[2].x);
		J(index, 4) = dy / (h / track[i].imgPt[2].y);
		J(index, 5) = dy / h;
		index++;

		//printf("%g %g %g\n", track[i].image - track[i].synth_image, dx / h, dy / h);
	}

	return J;
}

Matrix<double> generate_error_vector(vector<AlignedLOLAShot> & track, int num_points)
{
	Matrix<double> r(num_points, 1);
	int index = 0;
	for (unsigned int i = 0; i < track.size(); i++)
		if (track[i].image != -1 && track[i].synth_image != -1 && index < num_points)
		{
			r(index, 0) = track[i].image - track[i].synth_image;
			index++;
		}
	return r;
}

/**
 * B = transformation matrix
 * f_i(B) = synthetic image at track point i given B
 * r_i(B) = error at track point i given B
 * B_{s+1} = B_s + \delta
 * (J_f^T J_f) \delta = J_f^T r
 * Must compute J_f(B)
 * J_f(B)_{ij} = dF_i / dB_j
 **/
// assumes we have already called transform_track to set image and synth_image
Matrix3x3 gauss_newton_track(vector<AlignedLOLAShot> & track, ImageView<PixelGray<float> > cubImage,
		Matrix3x3 matrix)
{
	Matrix3x3 B = matrix;

	int num_points = 0;
	float err, last_err = compute_transform_error(track, &num_points);
	if (num_points <= 0) // no points
		return matrix;
	while (true)
	{
		Matrix<double> J = compute_jacobian(track, B, cubImage, num_points);
		Matrix<double> trans = transpose(J);
		Matrix<double> r = generate_error_vector(track, num_points);
		Matrix<double> U(6, 6), S(6, 6), VT(6, 6);
		Vector<double> s(6);
		svd(trans * J, U, s, VT);
		for (int i = 0; i < 6; i++)
		{
			if (s(i) > 10e-15) // don't divide by 0
				S(i, i) = 1.0 / s(i);
		}
		Matrix<double> pseudoinverse = transpose(VT) * S * transpose(U);
		Matrix<double> delta = pseudoinverse * trans * r;
		Matrix3x3 last_matrix = B;
		B(0, 0) += delta(0, 0);
		B(0, 1) += delta(1, 0);
		B(0, 2) += delta(2, 0);
		B(1, 0) += delta(3, 0);
		B(1, 1) += delta(4, 0);
		B(1, 2) += delta(5, 0);
		transform_track(track, B, cubImage);
		err = compute_transform_error(track);
		if (err == -1)
		{
			fprintf(stderr, "Matrix jumped too far, giving up.\n");
			return last_matrix;
		}
		//printf("Error: %g Old Error: %g\n", err, last_err);
		if (err >= last_err)
		{
			B = last_matrix;
			break;
		}
		last_err = err;
	}
	return B;
}


#define CHUNKSIZE   1


void SaveReportFile(vector<vector<LOLAShot> > &trackPts, vector<Vector<float, 6> >finalTransfArray,  
                      vector<float> errorArray, string reportFilename)
{
    FILE *fp = fopen(reportFilename.c_str(),"w");
    int numElements = errorArray.size();
    int bestResult = 0;
    float smallestError = 1000000000.0;

    for (int i = 0; i < numElements; i++){
       fprintf(fp,"index=%d d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f e=%f\n", 
	               i,
                       finalTransfArray[i](0), finalTransfArray[i](1), 
                       finalTransfArray[i](2), finalTransfArray[i](3),
                       finalTransfArray[i](4), finalTransfArray[i](5),
	               errorArray[i]);
        if ( (errorArray[i] < smallestError) && (errorArray[i] > 0) ){
          smallestError = errorArray[i]; 
          bestResult = i;
        }    
      
    }
    //print the best transform index and smallest error
    fprintf(fp, "bestTransformIndex = %d, smallest error = %f\n", bestResult, errorArray[bestResult]);
    //print the total number of tracks
    fprintf(fp, "numTotalTracks = %d\n", (int)trackPts.size());
    
    //print the number of features per track
    for (unsigned int ti = 0; ti < trackPts.size(); ti++){
      fprintf(fp,"trackID: %d ", ti);
      fprintf(fp,"numShots: %d ", (int)trackPts[ti].size());
      int numFeatures = 0;
      for (unsigned int si = 0; si < trackPts[ti].size(); si++){
	if ((trackPts[ti][si].featurePtRefl == 1.0) && (trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0)) { 
	  numFeatures++;  
	}
      }
      fprintf(fp,"numFeatures: %d\n", numFeatures);
   }


   fclose(fp);
}




/*
bool deriv_cached(string & input_name, string & cached_table){
  bool identified = false;
  string line;
  fstream myfile(cached_table.c_str());
  // does the name exist
  if(myfile.is_open()){
    while(! myfile.eof()){
      getline(myfile,line);
      cout << line << endl;
      if(line == input_name){
        identified = true;
        cout << "Fount it!" << endl; 
      }
    }
  }
  myfile.close();

  //if not add the name to the cached file
  if(!identified){
    // if we didn't find the file append the file to the document of saves
    ofstream file_change;
    file_change.open(cached_table.c_str(),ios::out | ios:: app);
    file_change << input_name << "\n";
    file_change.close();
  }


  return identified;
}
*/
void printRHS_LHS( Matrix<float,6,6> rhs,  Vector<float,6> lhs, int index, int iter)
{
  cout<<"index= "<<index<<", iter= "<<iter<<endl;
  printf("--------------------------------------\n");
  printf("RHS = \n");
  printf("[ %f %f %f | %f %f %f]\n",rhs(0,0), rhs(0,1), rhs(0,2), rhs(0,3), rhs(0,4), rhs(0,5) );   
  printf("[ %f %f %f | %f %f %f]\n", rhs(1,0), rhs(1,1), rhs(1,2), rhs(1,3), rhs(1,4),rhs(1,5) ); 
  printf("[ %f %f %f | %f %f %f]\n", rhs(2,0), rhs(2,1),rhs(2,2), rhs(2,3), rhs(2,4), rhs(2,5) );
  printf("[ %f %f %f | %f %f %f]\n", rhs(3,0), rhs(3,1), rhs(3,2), rhs(3,3), rhs(3,4), rhs(3,5) );
  printf("[ %f %f %f | %f %f %f]\n", rhs(4,0), rhs(4,1), rhs(4,2), rhs(4,3), rhs(4,4), rhs(4,5) );
  printf("[ %f %f %f | %f %f %f]\n", rhs(5,0), rhs(5,1), rhs(5,2), rhs(5,3), rhs(5,4), rhs(5,5) );
  printf("--------------------------------------\n\n");
  printf("LHS = [ %f, %f, %f, %f, %f, %f]\n", lhs[0], lhs[1], lhs[2], lhs[3], lhs[4], lhs[5]);
  printf("--------------------------------------\n\n");
}

void printIterationValues(  float error, int index, int iter, int numValidPts)
{
  cout<<"index= "<<index<<", iter= "<<iter<<", error="<<error<<", numValidPts="<<numValidPts<<endl;
}
/*
//This functions needs to be changed to generate multiple starting points for each 
//element of the affine transform 
void GenerateInitTransforms( vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings)
{
    initTransfArray.resize(settings.maxNumStarts);
    for (int i = 0; i < settings.maxNumStarts; i++){
      initTransfArray[i][0] = 1.0;
      initTransfArray[i][1] = 0.0;
      initTransfArray[i][2] = (i-settings.maxNumStarts/2)*2;
      initTransfArray[i][3] = 0.0;
      initTransfArray[i][4] = 1.0;
      initTransfArray[i][5] = 0.0;//(i-maxNumStarts/2)*25;
    }  
}
*/

//determines the matches between lidar tracks (trackPts) and images (cubFilename) around the lidar feature points.
//matchWindowHalfSize is half the size of the search window in image domain (pixels).
//numSamples is the number of LOLA points considered in matching (one point is unreliable due to image noise, 
//all points will violate the orthoprojection assumption).
vector<Vector4> FindMatches2D(vector<vector<AlignedLOLAShot> > &trackPts, string cubFilename, 
                              Vector2 matchWindowHalfSize, int numSamples, vector<float> &errorArray)
{
    vector<Vector4> matchArray;
    //vector<float>   errorArray;
   
    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
    double noDataValue = rsrc->nodata_read();

    //DiskImageView<PixelGray<float> > img( rsrc );
    DiskImageView<float>  img( rsrc );
    int width = img.cols();
    int height = img.rows();

    //Vector2 gain_bias = ComputeGainBiasFactor(trackPts);
    
    int matchWindowHalfWidth = matchWindowHalfSize(0);
    int matchWindowHalfHeight = matchWindowHalfSize(1);

    for (unsigned int ti = 0; ti < trackPts.size(); ti++){//for each track
      Vector2 gain_bias = ComputeGainBiasFactor(trackPts[ti]);
      for (unsigned int si = 0; si < trackPts[ti].size(); si++){//for each shot
	if ((trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1) && (trackPts[ti][si].featurePtLOLA == 1)){
          //valid track, reflectance and featurePt
	  
          float minDist = 10000000.0;
	  Vector2 bestMatch;
	  bestMatch(0) = 0;
	  bestMatch(1) = 0;
	  
	  int firstSample = si-numSamples/2;
	  if (firstSample < 0){
	    firstSample = 0;
	  }
	  unsigned int lastSample = si+numSamples/2;
	  if (lastSample > trackPts[ti].size()-1){
	    lastSample = trackPts[ti].size()-1;
	  }

	  
	  for (int k =  -matchWindowHalfHeight; k < matchWindowHalfHeight+1; k++){
	    for (int l = -matchWindowHalfWidth; l < matchWindowHalfWidth+1; l++){

              float dist = 0.0;   
              int numValidSamples = 0;

	      for (unsigned int i = firstSample; i< lastSample; i++){

		int x = (int)floor(trackPts[ti][i].imgPt[0].x); 
		int y = (int)floor(trackPts[ti][i].imgPt[0].y); 
		 
		if ((x+l > 0) && (y+k > 0) && (x+l < width) && (y+k < height)){
		  if ((img(x+l,y+k)!=noDataValue) && (gain_bias(0)!= 0) && (gain_bias(1)!=0)){
		
		    dist = dist + fabs(img(x+l,y+k) - gain_bias(1) - gain_bias(0)*trackPts[ti][i].reflectance);
		    //dist = dist + fabs(img(x+l,y+k) - trackPts[ti][i].reflectance);
                    numValidSamples++;
		  }
		}

	      }//i

              //at least half the samples are valid
              if (numValidSamples > numSamples/2){ 
                 dist = dist/numValidSamples;
		 if (dist < minDist){
		   minDist = dist;
		   bestMatch(0) = l;
		   bestMatch(1) = k;
		   //cout<<"minDist="<<minDist<<", l="<<l<<", k="<<k<<endl;
		 }      
              }


	    }//k	
	  }//l
	    
          Vector4 feature_match;
	  feature_match(0) = trackPts[ti][si].imgPt[0].x;
	  feature_match(1) = trackPts[ti][si].imgPt[0].y;
	  feature_match(2) = feature_match(0) + bestMatch(0);
	  feature_match(3) = feature_match(1) + bestMatch(1);
	  matchArray.push_back(feature_match);
          errorArray.push_back(minDist);
          
        

	}
      }
    }

    return matchArray;
}

// trackPts is vector of tracks
/*vector<Vector4> FindMatches2D(vector<vector<LOLAShot> > &trackPts, string cubFilename, Vector2 matchWindowHalfSize)
{
    vector<Vector4> matchArray;
   
    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
    double noDataValue = rsrc->nodata_read();

    DiskImageView<PixelGray<float> > img( rsrc );
    int width = img.cols();
    int height = img.rows();

    //Vector2 gain_bias = ComputeGainBiasFactor(trackPts);
    
    int matchWindowHalfWidth = matchWindowHalfSize(0);
    int matchWindowHalfHeight = matchWindowHalfSize(1);

    for (unsigned int ti = 0; ti < trackPts.size(); ti++){//for each track
      for (unsigned int si = 0; si < trackPts[ti].size(); si++){//for each shot
	if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) &&(trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
	  
	  int x = (int)floor(trackPts[ti][si].imgPt[0].x); 
	  int y = (int)floor(trackPts[ti][si].imgPt[0].y); 
	  
	  float minDist = 10000000.0;
	  Vector2 bestMatch;
	  bestMatch(0) = x;
	  bestMatch(1) = y;
	
	  for (int k =  y - matchWindowHalfHeight; k < y + matchWindowHalfHeight+1; k++){
	    for (int l = x - matchWindowHalfWidth; l < x + matchWindowHalfWidth+1; l++){
	      if ((l > 0) && (k > 0) && (l < width) && (k < height)){
		if (img(l,k)!=noDataValue){
		  //float dist = fabs(img(l,k) - gain_bias(1) - gain_bias(0)*trackPts[ti][si].reflectance);
		  float dist = fabs(img(l,k) - trackPts[ti][si].reflectance);
		  if (dist < minDist){
		    minDist = dist;
		    bestMatch(0) = l;
		    bestMatch(1) = k;
		  }     
		}
	      }
	    }
	  }	
  
          Vector4 feature_match;
          feature_match(0)=trackPts[ti][si].imgPt[0].x;
          feature_match(1)=trackPts[ti][si].imgPt[0].y;
          feature_match(2)=bestMatch(0);
          feature_match(3)=bestMatch(1);
          matchArray.push_back(feature_match);
        
	}
      }
    }

    return matchArray;
}


void EstimateMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename, Vector2 matchWindowHalfSize, 
                                   vector<Vector<float, 6> > &initTransfArray, vector<float> &matchingErrorArray)
{
 

 cout<<"Estimate matching params from cub..."<<endl;
 

 boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );

 DiskImageView<PixelGray<float> > img( rsrc );

 //Vector2 gain_bias = ComputeGainBiasFactor(trackPts);

 cout<<"Computing LOLA centroid ..."<<endl;
 float i_C = 0.0;
 float j_C = 0.0;
 int numPts = 0;
 
 //compute the  i_C and  j_C new way - START
 for (unsigned int ti = 0; ti < trackPts.size(); ti++){
   for (unsigned int si = 0; si < trackPts[ti].size(); si++){
     if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
       //if ((float)img(trackPts[ti][si].imgPt[0].x, trackPts[ti][si].imgPt[0].y)>minmax(0)){
	 i_C = i_C + trackPts[ti][si].imgPt[0].x;
	 j_C = j_C + trackPts[ti][si].imgPt[0].y;
	 numPts++;
	 //}
     }
   }
 }
 
 i_C = i_C/numPts;
 j_C = j_C/numPts;

 float sum_x2 = 0.0;
 float sum_xy = 0.0;
 float sum_x  = 0.0;
 float sum_y  = 0.0;
 float sum_y2 = 0.0;
 float sum_1  = 0.0;

 float sum_kx = 0.0; 
 float sum_ky = 0.0;
 float sum_k  = 0.0;
 float sum_lx = 0.0;
 float sum_ly = 0.0;
 float sum_l  = 0.0;

 vector<float> errorArray;
 vector<Vector4> matchArray = FindMatches2D(trackPts, cubFilename, matchWindowHalfSize, 20, errorArray);

 cout<<"num_matches="<<matchArray.size()<<endl;


 for (unsigned int i = 0; i < matchArray.size(); i++){

       Vector2 bestMatch;
       float x, y;

       x = matchArray[i](0)-i_C;
       y = matchArray[i](1)-j_C;
       bestMatch(0) = matchArray[i](2)-i_C;
       bestMatch(1) = matchArray[i](3)-j_C;
       cout<<"error="<<errorArray[i]<<endl;       

       sum_x2 = sum_x2 + x*x;
       sum_xy = sum_xy + x*y;
       sum_x  = sum_x  + x;
       sum_y  = sum_y  + y;
       sum_y2 = sum_y2 + y*y;
       sum_1  = sum_1 + 1;

       sum_lx = sum_lx + bestMatch(0)*x;
       sum_ly = sum_ly + bestMatch(0)*y;
       sum_l  = sum_l  + bestMatch(0);
       
       sum_kx = sum_kx + bestMatch(1)*x; 
       sum_ky = sum_ky + bestMatch(1)*y;
       sum_k  = sum_k  + bestMatch(1);
 }

 //compute the affine transform
 cout<<"numsamples = "<<sum_1<<endl;
 Matrix<float,3,3> rhs;
 Vector<float,3> lhs;
 initTransfArray.resize(1);

 rhs(0,0) = sum_x2/sum_1;
 rhs(0,1) = sum_xy/sum_1;
 rhs(0,2) = sum_x/sum_1;

 rhs(1,0) = sum_xy/sum_1;
 rhs(1,1) = sum_y2/sum_1;
 rhs(1,2) = sum_y/sum_1;

 rhs(2,0) = sum_x/sum_1;
 rhs(2,1) = sum_y/sum_1;
 rhs(2,2) = sum_1/sum_1;

 lhs(0) = sum_lx/sum_1;
 lhs(1) = sum_ly/sum_1;
 lhs(2) = sum_l/sum_1;
 cout<<"rhs"<<rhs<<endl;
 cout<<"lhs"<<lhs<<endl;
 
 try {
        solve_symmetric_nocopy(rhs,lhs);
        initTransfArray[0][0] = lhs(0); 
	initTransfArray[0][1] = lhs(1);
	initTransfArray[0][2] = lhs(2);
      } catch (ArgumentErr &) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
      }

 lhs(0) = sum_kx/sum_1;
 lhs(1) = sum_ky/sum_1;
 lhs(2) = sum_k/sum_1;

 rhs(0,0) = sum_x2/sum_1;
 rhs(0,1) = sum_xy/sum_1;
 rhs(0,2) = sum_x/sum_1;

 rhs(1,0) = sum_xy/sum_1;
 rhs(1,1) = sum_y2/sum_1;
 rhs(1,2) = sum_y/sum_1;

 rhs(2,0) = sum_x/sum_1;
 rhs(2,1) = sum_y/sum_1;
 rhs(2,2) = sum_1/sum_1;

 try {
   solve_symmetric_nocopy(rhs,lhs);
   initTransfArray[0][3] = lhs(0); 
   initTransfArray[0][4] = lhs(1);
   initTransfArray[0][5] = lhs(2);
 } catch (ArgumentErr &e) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
 }



 cout<<initTransfArray[0]<<endl;
 matchingErrorArray.resize(1);
 matchingErrorArray[0] = 1000000.0;
 matchArray.clear();

};*/

void EstimateAffineTransform(vector<Vector<float, 6> > &initTransfArray, vector<vector<LOLAShot> > &trackPts, 
                             vector<Vector4> matchArray, vector<float> errorArray)
{
 

 cout<<"Estimate affine transform params ..."<<endl;
 
 /*
 boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
 double noDataValue = rsrc->nodata_read();

 DiskImageView<PixelGray<float> > img( rsrc );
 int width = img.cols();
 int height = img.rows();

 Vector2 gain_bias = ComputeGainBiasFactor(trackPts);
 */
 cout<<"Computing LOLA centroid ..."<<endl;
 float i_C = 0.0;
 float j_C = 0.0;
 int numPts = 0;
 
 //compute the  i_C and  j_C new way - START
 for (unsigned int ti = 0; ti < trackPts.size(); ti++){
   for (unsigned int si = 0; si < trackPts[ti].size(); si++){
     if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
       //if ((float)img(trackPts[ti][si].imgPt[0].x, trackPts[ti][si].imgPt[0].y)>minmax(0)){
	 i_C = i_C + trackPts[ti][si].imgPt[0].x;
	 j_C = j_C + trackPts[ti][si].imgPt[0].y;
	 numPts++;
	 //}
     }
   }
 }
 
 i_C = i_C/numPts;
 j_C = j_C/numPts;


 float sum_x2 = 0.0;
 float sum_xy = 0.0;
 float sum_x  = 0.0;
 float sum_y  = 0.0;
 float sum_y2 = 0.0;
 float sum_1  = 0.0;

 float sum_kx = 0.0; 
 float sum_ky = 0.0;
 float sum_k  = 0.0;
 float sum_lx = 0.0;
 float sum_ly = 0.0;
 float sum_l  = 0.0;
 /*
 int matchWindowHalfWidth = matchWindowHalfSize(0);
 int matchWindowHalfHeight = matchWindowHalfSize(1);
 vector<float> errorArray;
 vector<Vector4> matchArray = FindMatches2D(trackPts, cubFilename, matchWindowHalfSize, 20, errorArray);
 */
 cout<<"num_matches="<<matchArray.size()<<endl;


 for (unsigned int i = 0; i < matchArray.size(); i++){

       Vector2 bestMatch;
       float x, y;

       x = matchArray[i](0)-i_C;
       y = matchArray[i](1)-j_C;
       bestMatch(0) = matchArray[i](2)-i_C;
       bestMatch(1) = matchArray[i](3)-j_C;
       //cout<<"error="<<errorArray[i]<<endl;       

       sum_x2 = sum_x2 + x*x;
       sum_xy = sum_xy + x*y;
       sum_x  = sum_x  + x;
       sum_y  = sum_y  + y;
       sum_y2 = sum_y2 + y*y;
       sum_1  = sum_1 + 1;

       sum_lx = sum_lx + bestMatch(0)*x;
       sum_ly = sum_ly + bestMatch(0)*y;
       sum_l  = sum_l  + bestMatch(0);
       
       sum_kx = sum_kx + bestMatch(1)*x; 
       sum_ky = sum_ky + bestMatch(1)*y;
       sum_k  = sum_k  + bestMatch(1);
 }

 //compute the affine transform
 cout<<"numsamples = "<<sum_1<<endl;
 Matrix<float,3,3> rhs;
 Vector<float,3> lhs;
 initTransfArray.resize(1);

 rhs(0,0) = sum_x2/sum_1;
 rhs(0,1) = sum_xy/sum_1;
 rhs(0,2) = sum_x/sum_1;

 rhs(1,0) = sum_xy/sum_1;
 rhs(1,1) = sum_y2/sum_1;
 rhs(1,2) = sum_y/sum_1;

 rhs(2,0) = sum_x/sum_1;
 rhs(2,1) = sum_y/sum_1;
 rhs(2,2) = sum_1/sum_1;

 lhs(0) = sum_lx/sum_1;
 lhs(1) = sum_ly/sum_1;
 lhs(2) = sum_l/sum_1;
 
 try {
        solve_symmetric_nocopy(rhs,lhs);
        initTransfArray[0][0] = lhs(0); 
	initTransfArray[0][1] = lhs(1);
	initTransfArray[0][2] = lhs(2);
      } catch (ArgumentErr &/*e*/) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
      }

 lhs(0) = sum_kx/sum_1;
 lhs(1) = sum_ky/sum_1;
 lhs(2) = sum_k/sum_1;

 rhs(0,0) = sum_x2/sum_1;
 rhs(0,1) = sum_xy/sum_1;
 rhs(0,2) = sum_x/sum_1;

 rhs(1,0) = sum_xy/sum_1;
 rhs(1,1) = sum_y2/sum_1;
 rhs(1,2) = sum_y/sum_1;

 rhs(2,0) = sum_x/sum_1;
 rhs(2,1) = sum_y/sum_1;
 rhs(2,2) = sum_1/sum_1;

 try {
   solve_symmetric_nocopy(rhs,lhs);
   initTransfArray[0][3] = lhs(0); 
   initTransfArray[0][4] = lhs(1);
   initTransfArray[0][5] = lhs(2);
 } catch (ArgumentErr &/*e*/) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
 }



 cout<<initTransfArray[0]<<endl;
 //matchingErrorArray.resize(1);
 //matchingErrorArray[0] = 1000000.0;
 //matchArray.clear();

};

//This functions needs to be changed to generate multiple starting points for each 
//element of the affine transform 
void GenerateInitTransforms(vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings)
{

  vector<vector<float> > params;
  params.resize(6);
  
  params[0].resize(3);
  for (unsigned int i = 0; i < params[0].size(); i++){
    params[0][i] = (float)(0.995+i*0.005);
    //params[0][i] = (float)(0.975+i*0.025);
  }
  
  /*
  params[0].resize(1);
  params[0][0] = 1;  
  */
  
  params[1].resize(3);
  for (unsigned int i = 0; i < params[1].size(); i++){
    params[1][i] = (float)(-0.005+i*0.005);
    //params[1][i] = (float)(-0.025+i*0.025);
  }
  /*
  params[1].resize(1);
  params[1][0] = 0;
  */
  params[2].resize(9);
  for (unsigned int i = 0; i < params[2].size(); i++){
    params[2][i] = (float)(-8.0+2*i);
  }
 
  
  params[3].resize(3);
  for (unsigned int i = 0; i < params[3].size(); i++){
    params[3][i] = (float)(-0.005+i*0.005);
    //params[3][i] = (float)(-0.025+i*0.025);
  }
  /*
  params[3].resize(1);
  params[3][0] = 0;  
  */
  
  params[4].resize(3);
  for (unsigned int i = 0; i < params[4].size(); i++){
    params[4][i] = (float)(0.995+i*0.005);
    //params[4][i] = (float)(0.975+i*0.025);
  }
  
  /*
  params[4].resize(1);
  params[4][0] = 1;
  */
  params[5].resize(9);
  for (unsigned int i = 0; i < params[5].size(); i++){
    params[5][i] = (float)(-8.0+2*i);
  }
 

  int maxNumStarts = 1;
  for( int i = 0; i < 6; i++){
    maxNumStarts = maxNumStarts*params[i].size();
  } 
  initTransfArray.resize(maxNumStarts);
  
  int index = 0;
  for (unsigned int i0 = 0; i0 < params[0].size(); i0++){
    for (unsigned int i1 = 0; i1 < params[1].size(); i1++){ 
      for (unsigned int i2 = 0; i2 < params[2].size(); i2++){
        for (unsigned int i3 = 0; i3 < params[3].size(); i3++){
          for (unsigned int i4 = 0; i4 < params[4].size(); i4++){
            for (unsigned int i5 = 0; i5 < params[5].size(); i5++){
	      initTransfArray[index][0] = params[0][i0];
	      initTransfArray[index][1] = params[1][i1];
	      initTransfArray[index][2] = params[2][i2];
	      initTransfArray[index][3] = params[3][i3];
	      initTransfArray[index][4] = params[4][i4];
	      initTransfArray[index][5] = params[5][i5];
              index++;
	    }
	  }
	}
      }
    }
  }
  /*
  for (int i = 0; i < maxNumStarts; i++){
    printf("%d:", i);
    for (int j = 0; j < 6; j++){
      printf("%f ", initTransfArray[i][j]);
    }
    printf("\n");
  }
  */
}

void GetBestTransform(vector<Vector<float, 6> > &finalTransfArray,  vector<float> &finalMatchingErrorArray, 
                      Vector<float, 6> &optimalTransf, float &optimalError)
{
    //this will be made a special function to make sure it is flexible wrt to the name of the parameters.
    int bestResult = 0;
    float smallestError = std::numeric_limits<float>::max();
    for (unsigned int index = 0; index < finalTransfArray.size(); index++){
      printf("refined %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
	     index, finalMatchingErrorArray[index], 
	     finalTransfArray[index](0), finalTransfArray[index](1), 
	     finalTransfArray[index](2), finalTransfArray[index](3),
	     finalTransfArray[index](4), finalTransfArray[index](5));
      if  ( (finalMatchingErrorArray[index] < smallestError) && (finalMatchingErrorArray[index] > 0) ){
	smallestError = finalMatchingErrorArray[index]; 
	bestResult = index;
      }    
    } 
    cout<<"bestResult= "<<bestResult<<endl;

    //copy the best transform to optimalTransfArray
    optimalTransf = finalTransfArray[bestResult];
    optimalError  = finalMatchingErrorArray[bestResult];
}

//determines the best finalTransfArray from all initTransfArrays
//it assumes that each image is transformed by one affine transform
//in the future we can investigate the use of one affine transform per track. 
/*void InitMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename,  
			       CoregistrationParams coregistrationParams,  
			       vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> > &finalTransfArray, 
			       vector<float> &matchingErrorArray)
{
    unsigned int ti, si, li;
    float jA, iA;
    float I_e_val;
    Vector2 minmax = ComputeMinMaxValuesFromCub(cubFilename);
    cout<<"minmax="<<minmax<<endl;

   
    vector<float> errorArray;
    errorArray.resize(initTransfArray.size());
    vector<int> numValidPts;
    numValidPts.resize(initTransfArray.size());    

    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  
    DiskImageView<PixelGray<float> > isis_view( rsrc );
    int width = isis_view.cols();
    int height = isis_view.rows();
  
    camera::IsisCameraModel model(cubFilename);

    InterpolationView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, BilinearInterpolation> interpImg
    = interpolate(isis_view, BilinearInterpolation(), ConstantEdgeExtension());
  
    unsigned int index;
    int initTransformIndex = 0;

    //compute the scale factor between LOLA reflectance and image.
    cout<<"Computing the gain and bias ..."<<endl;
    Vector2 gain_bias = ComputeGainBiasFactor(trackPts);
    cout<<"gain_bias="<<gain_bias<<endl;   

    cout<<"Computing LOLA centroid ..."<<endl;
    float i_C = 0.0;
    float j_C = 0.0;
    int numPts = 0;

    //compute the  i_C and  j_C new way - START
    for (unsigned int ti = 0; ti < trackPts.size(); ti++){
      for (unsigned int si = 0; si < trackPts[ti].size(); si++){
	if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
	  for (unsigned int li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot
	    if (trackPts[ti][si].LOLAPt[li].s == 1){//center point of a valid shot
	      if ((float)interpImg(trackPts[ti][si].imgPt[li].x, trackPts[ti][si].imgPt[li].y)>minmax(0)){
		i_C = i_C + trackPts[ti][si].imgPt[li].x;
		j_C = j_C + trackPts[ti][si].imgPt[li].y;
		numPts++;
	      }
	    }
	  }
	}
      }
    }
    
    i_C = i_C/numPts;
    j_C = j_C/numPts;

    cout <<"done."<<endl;

    for (index = 0; index < initTransfArray.size(); index++){
      
      errorArray[index] = 0.0;
      numValidPts[index] = 0;
      if ((initTransfArray[index][2]==0.0) && (initTransfArray[index][5]==0.0)){
	initTransformIndex = index;
      }

      for (ti = 0; ti < trackPts.size(); ti++){//for each track

	for (si = 0; si < trackPts[ti].size(); si++){//for each shot
          
          //if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0)){
          if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) &&(trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
           
            //float weight = trackPts[ti][si].weightRefl;
        
            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 1){//center point of a valid shot
    
                 
		iA = initTransfArray[index][0]*(trackPts[ti][si].imgPt[li].x-i_C) 
		   + initTransfArray[index][1]*(trackPts[ti][si].imgPt[li].y-j_C) 
                   + initTransfArray[index][2]+i_C;
		jA = initTransfArray[index][3]*(trackPts[ti][si].imgPt[li].x-i_C) 
		   + initTransfArray[index][4]*(trackPts[ti][si].imgPt[li].y-j_C) 
                   + initTransfArray[index][5]+j_C;

                            
                //check (iA,jA) are inside the image!
		if ( ( iA >= 1) && ( iA < height-1) && ( jA >= 1) && ( jA < width-1)){   
                  if ((float)interpImg(jA, iA)> minmax(0)){//valid value  
                    I_e_val = (float)interpImg(jA,iA) - gain_bias(1) - gain_bias(0)*trackPts[ti][si].reflectance;
		    errorArray[index] += abs(I_e_val);
		    numValidPts[index]++;
		  }
                  else{
                    I_e_val= 0;
                  }
                  //cout<<"I_e_val="<<I_e_val<<endl;
		  //errorArray[index] += abs(I_e_val);
                  //numValidPts[index]++;
		}
		
	      }
	    }
	  }
	}
      }
    }

    //determine the best image transform for which the matching error is minimized 
    float minError = errorArray[0]/numValidPts[0];
    int bestIndex = 0;
    finalTransfArray.resize(1);
    matchingErrorArray.resize(1);
    finalTransfArray[0] = initTransfArray[0];

    for (index = 0; index < initTransfArray.size(); index++){
      errorArray[index] = errorArray[index]/numValidPts[index];  
      if (errorArray[index] < minError){
	 minError = errorArray[index];
         finalTransfArray[0] = initTransfArray[index];
         bestIndex = index;
      }  
    }
    cout<<"initError: "<<errorArray[initTransformIndex]<<endl;
    cout <<"finalTransform: "<<finalTransfArray[0]<<endl;
    cout<<"minError= "<<minError<<endl;
    cout<<"bestIndex="<<bestIndex<<endl; 

    matchingErrorArray[0] = minError;
    errorArray.clear();
    numValidPts.clear();  
    
}*/

//image to lidar coregistration
/*void UpdateMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename,  
			         int numMaxIter, 
			         vector<Vector<float, 6> >initTransfArray, vector<float> &initErrorArray, 
                                 vector<Vector<float, 6> >&finalTransfArray, vector<float> &errorArray, Vector2 &centroid)
{

  printf("LIMA matching single processor\n");

  Vector2 minmax = ComputeMinMaxValuesFromCub(cubFilename);

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodata_value = rsrc->nodata_read();
  cout<<"nodata value="<<nodata_value<<endl;
   
  DiskImageView<PixelGray<float> > isis_view( rsrc );
  int width = isis_view.cols();
  int height = isis_view.rows();
   
  cout <<"Interpolating the image ...";
  InterpolationView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, BilinearInterpolation> interpImg
    = interpolate(isis_view, BilinearInterpolation(), ConstantEdgeExtension());
  cout<<"done."<<endl;

  cout << "Computing the gain and bias...";
  Vector2 gain_bias = ComputeGainBiasFactor(trackPts);
  cout<<"done."<<endl;

  cout<<"width="<<width<<" height="<<height<<" gain="<<gain_bias(0)<<", bias="<<gain_bias(1)<<endl;

  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  //std::string temp = sufix_from_filename(cubFilename);
  std::string temp = GetFilenameNoPath(cubFilename);

  std::string xDerivFilename = "../aux/" + GetFilenameNoExt(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../aux/" + GetFilenameNoExt(temp) + "_y_deriv.tif";

  cout<<xDerivFilename<<endl;
  cout<<yDerivFilename<<endl;

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(isis_view),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float>x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(isis_view),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;

  cout <<"Interpolating the image derivatives...";
  InterpolationView<EdgeExtensionView<DiskImageView<float>, ConstantEdgeExtension>, BilinearInterpolation> interp_x_deriv
  = interpolate(x_deriv, BilinearInterpolation(), ConstantEdgeExtension());
  InterpolationView<EdgeExtensionView<DiskImageView<float>, ConstantEdgeExtension>, BilinearInterpolation> interp_y_deriv
    = interpolate(y_deriv, BilinearInterpolation(), ConstantEdgeExtension());
  cout<<"done."<<endl;
  
  cout<<"Computing LOLA centroid ..."<<endl;
  float i_C = 0.0;
  float j_C = 0.0;
  int numPts = 0;

  //compute the  i_C and  j_C new way - START
  for (unsigned int ti = 0; ti < trackPts.size(); ti++){
   for (unsigned int si = 0; si < trackPts[ti].size(); si++){
     if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
       for (unsigned int li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot
	 if (trackPts[ti][si].LOLAPt[li].s == 1){//center point of a valid shot
	   if ((float)interpImg(trackPts[ti][si].imgPt[li].x, trackPts[ti][si].imgPt[li].y)>minmax(0)){
               i_C = i_C+ trackPts[ti][si].imgPt[li].x;
               j_C = j_C+ trackPts[ti][si].imgPt[li].y;
               numPts++;
	    }
	  }
        }
      }
    }
  }

  i_C = i_C/numPts;
  j_C = j_C/numPts;
  centroid(0) = i_C;
  centroid(1) = j_C;
  cout <<"done."<<endl;

  int row_max, col_max;
  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
    
  errorArray.resize(initTransfArray.size());
  finalTransfArray.resize(initTransfArray.size());
  vector<int> numValidPts;
  numValidPts.resize(initTransfArray.size());
  
  //copy initTransfArray into finalTransfArray.
  for (unsigned int index = 0; index < initTransfArray.size(); index++){
    for (unsigned int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
    cout <<"InitTransform: "<<initTransfArray[index]<<endl;
  }

  for (unsigned int index = 0; index < initTransfArray.size(); index++){
  
    unsigned int ti,si,li;//indices for tracks, shots and lola points respectively
    float iA, jA;
    int ii, jj;
    int iter;
    float I_x_sqr, I_x_I_y, I_y_sqr; 
    float I_y_val, I_x_val;
    float I_e_val ;
    Matrix<float,6,6> rhs;
    Vector<float,6> lhs;

    iA = 0;
    jA = 0;
    iter = 0;

    while( iter <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs(i_RHS,j_RHS) = 0.0;
        }
        lhs(i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;
      numValidPts[index] = 0;
      
      for (ti = 0; ti < trackPts.size(); ti++){

        for (si = 0; si < trackPts[ti].size(); si++){

	  if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
           
            float weight = trackPts[ti][si].weightRefl;
            //cout<<"weight="<<weight<<endl;        

            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 1){//center point of a valid shot
    
		iA = finalTransfArray[index][0]*(trackPts[ti][si].imgPt[li].x-i_C)+ 
                     finalTransfArray[index][1]*(trackPts[ti][si].imgPt[li].y-j_C)+
                     finalTransfArray[index][2]+i_C;
		jA = finalTransfArray[index][3]*(trackPts[ti][si].imgPt[li].x-i_C)+ 
                     finalTransfArray[index][4]*(trackPts[ti][si].imgPt[li].y-j_C)+ 
                     finalTransfArray[index][5]+j_C;
                
                // check (iA,jA) are inside the image!
		if ( ( iA >= 1) && ( iA < row_max-1) && ( jA >= 1) && ( jA < col_max-1)){ 
	
                  if (((float)interpImg(jA,iA)>minmax(0)) && (fabs(interp_x_deriv(jA, iA))< 10000) && (fabs(interp_y_deriv(jA, iA))< 10000)){
 
		    // calculate ii & jj relative to the image center
		    
                    ii = iA - i_C;
		    jj = jA - j_C;
		    
		    I_e_val = (float)interpImg(jA,iA) - gain_bias(1) - gain_bias(0)*trackPts[ti][si].reflectance;
                     
                    //cout<<I_e_val<<endl;
                    //start here
                    
		    I_x_val = interp_x_deriv(jA,iA)*weight; 
		    I_y_val = interp_y_deriv(jA,iA)*weight; 
              
		    I_x_I_y = I_x_val*I_y_val;        
		    I_x_sqr = I_x_val*I_x_val;
		    I_y_sqr = I_y_val*I_y_val;

		    // Left hand side
		    lhs(0) += ii * I_x_val * I_e_val;
		    lhs(1) += jj * I_x_val * I_e_val;
		    lhs(2) +=      I_x_val * I_e_val;
		    lhs(3) += ii * I_y_val * I_e_val;
		    lhs(4) += jj * I_y_val * I_e_val;
		    lhs(5) +=      I_y_val * I_e_val;
		  
		    // Right Hand Side UL
		    rhs(0,0) += ii*ii * I_x_sqr;
		    rhs(0,1) += ii*jj * I_x_sqr;
		    rhs(0,2) += ii    * I_x_sqr;
		    rhs(1,1) += jj*jj * I_x_sqr;
		    rhs(1,2) += jj    * I_x_sqr;
		    rhs(2,2) +=         I_x_sqr;

		    // Right Hand Side UR
		    rhs(0,3) += ii*ii * I_x_I_y;
		    rhs(0,4) += ii*jj * I_x_I_y;
		    rhs(0,5) += ii    * I_x_I_y;
		    rhs(1,4) += jj*jj * I_x_I_y;
		    rhs(1,5) += jj    * I_x_I_y;
		    rhs(2,5) +=         I_x_I_y;

		    // Right Hand Side LR
		    rhs(3,3) += ii*ii * I_y_sqr;
		    rhs(3,4) += ii*jj * I_y_sqr;
		    rhs(3,5) += ii    * I_y_sqr;
		    rhs(4,4) += jj*jj * I_y_sqr;
		    rhs(4,5) += jj    * I_y_sqr;
		    rhs(5,5) +=         I_y_sqr;
		  //end here
                    
                    errorArray[index] += abs(I_e_val);
                    numValidPts[index]++;
		  }
                  else{
                    I_e_val= 0;
                  }

		  //errorArray[index] += abs(I_e_val);
                  //numValidPts[index]++;

		}
	      }           
            }// end of if statement: inside image
          }// end of if statement: valid reflectance  
        }// end of for loop over i
      }//end of loop over k

      // Fill in symmetric entries
      rhs(1,0) = rhs(0,1);
      rhs(2,0) = rhs(0,2);
      rhs(2,1) = rhs(1,2);
      rhs(1,3) = rhs(0,4);
      rhs(2,3) = rhs(0,5);
      rhs(2,4) = rhs(1,5);
      rhs(3,0) = rhs(0,3);
      rhs(3,1) = rhs(1,3);
      rhs(3,2) = rhs(2,3);
      rhs(4,0) = rhs(0,4);
      rhs(4,1) = rhs(1,4);
      rhs(4,2) = rhs(2,4);
      rhs(4,3) = rhs(3,4);
      rhs(5,0) = rhs(0,5);
      rhs(5,1) = rhs(1,5);
      rhs(5,2) = rhs(2,5);
      rhs(5,3) = rhs(3,5);
      rhs(5,4) = rhs(4,5);

      //printRHS_LHS(rhs, index, iter);
     
      try {
        solve_symmetric_nocopy(rhs,lhs);
      } catch (ArgumentErr &) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
      }

      //cout<<"numValidPts="<<numValidPts[index]<<endl; 
      errorArray[index] = errorArray[index]/numValidPts[index];
      printIterationValues(errorArray[index], index, iter, numValidPts[index]);

      finalTransfArray[index] += lhs;
  
      iter ++;
 
    }
   
    if (errorArray[index] > initErrorArray[index]){
       finalTransfArray[index] = initTransfArray[index];
       errorArray[index]  = initErrorArray[index];
    }
  }//index loop ends here

} */



#if 0
//determines the best finalTransfArray from all initTransfArrays
//it assumes that each image is transformed by one affine transform
//in the future we can investigate the use of one affine transform per track. 
void InitMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
			ModelParams modelParams, CoregistrationParams coregistrationParams,  
			vector<Vector<float, 6> >initTransfArray, Vector<float, 6> &finalTransf, 
			float *matchingError)
{
    unsigned int ti, si, li;
    int jA, iA;
    float I_e_val;
  
    vector<float>errorArray;
    errorArray.resize(initTransfArray.size());
    DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
    GeoReference DRGGeo;
    read_georeference(DRGGeo, DRGFilename);

    cout <<"Interpolating the image ...";
    ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
									 ConstantEdgeExtension()),
							     BilinearInterpolation());
    cout<<"done."<<endl;

    unsigned int index;
    float scaleFactor;
    int row_max = DRG.rows();
    int col_max = DRG.cols();
  
    //compute the scale factor between LOLA reflectance and image.
    scaleFactor = ComputeScaleFactor(trackPts);
    cout<<"done computing the scaling factor"<<endl;

    for (index = 0; index < initTransfArray.size(); index++){
      
      errorArray[index] = 0.0;
      
      for (ti = 0; ti < trackPts.size(); ti++){//for each track

	for (si = 0; si < trackPts[ti].size(); si++){//for each shot
          
          //cout<<ti<<" "<<si<<endl;

          if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != /*0*/-1)){
           
            // float weight = trackPts[ti][si].weightRefl;
        
            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 3){//center point of a valid shot
    
                 
		iA = (int) floor(initTransfArray[index][0]*trackPts[ti][si].imgPt[li].x 
                               + initTransfArray[index][1]*trackPts[ti][si].imgPt[li].y 
                               + initTransfArray[index][2]);
		jA = (int) floor(initTransfArray[index][3]*trackPts[ti][si].imgPt[li].x 
                               + initTransfArray[index][4]*trackPts[ti][si].imgPt[li].y 
                               + initTransfArray[index][5]);
                
                //check (iA,jA) are inside the image!
		if ( ( iA >= 0) && ( iA < row_max) && ( jA >= 0) && ( jA < col_max)){ 
		  if (is_valid(interpDRG(jA, iA)) ){
                    //calculate ii & jj relative to the image center
		    I_e_val = interpDRG(jA,iA) - scaleFactor*trackPts[ti][si].reflectance;
		  }
                  else{
                    I_e_val= 0;
                  }

		  errorArray[index] += abs(I_e_val);
		}
		
	      }
	    }
	  }
	}
      }
    }

    //determine the best image transform for which the matching error is minimized 
    float minError = errorArray[0];
    finalTransf = initTransfArray[0];
    for (index = 0; index < initTransfArray.size(); index++){
      if (errorArray[index] < minError){
	minError = errorArray[index];
        finalTransf = initTransfArray[index];
      }  
    } 
    *matchingError = minError;
}
#endif

#if 0
//image to lidar coregistration
void UpdateMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
			  ModelParams modelParams,  int numMaxIter, 
			  vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
			  vector<float> &errorArray )
{

  printf("LIMA matching single processor\n");

  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing the scale factor...";
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);
  cout<<"done."<<endl;

  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DRGFilename);

  std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_y_deriv.tif";

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  

  //float xx,yy;
  int row_max, col_max;

  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
   
  cout << "image size: rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C = row_max/2;
  int j_C = col_max/2;

  //working with the entire CVS file significant errors where seen at index 126 - set index to this by hand to debugg this test case
 /*
  1. we see huge problems from the get go - massive lhs numbers, where does these come from in comparison
  2. need to explore & instrument the code further
  3. we have some stable/unstable pairs to explore (106 stable, 109 not) & ( 112 stable, 115 not)
 
  A couple of thoughts - we rezero lhs & rhs every iteration this should not be the problem.  Could the problem simple be a scaling issue where the size of the error terms on the right & left hand side - caused by the great number of points - causes numerical stability issues for solving the matrix?
*/

  //copy initTransfArray into finalTransfArray.
  for (unsigned int index = 0; index < initTransfArray.size(); index++){
    for (unsigned int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
  }

  for (unsigned int index = 0; index < initTransfArray.size(); index++){

    cout << "index = "<< index << endl;
  
  
    unsigned int ti,si,li;//indices for tracks, shots and lola points respectively
    int iA, jA;
    int ii, jj;
    int iter;
    float I_x_sqr, I_x_I_y, I_y_sqr; 
    float I_y_val, I_x_val;
    float I_e_val ;
    Matrix<float,6,6> rhs;
    Vector<float,6> lhs;
 
    iA = 0;
    jA = 0;
    iter = 0;

    while( iter <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs(i_RHS,j_RHS) = 0.0;
        }
        lhs(i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;
      
      for (ti = 0; ti < trackPts.size(); ti++){

        for (si = 0; si < trackPts[ti].size(); si++){

          //if ((trackPts[ti][si].valid ==1) && (trackPts[ti][si].reflectance !=0)){
          if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) &&(trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
           
            float weight = trackPts[ti][si].weightRefl;
        
            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 3){//center point of a valid shot
    
		iA = (int) floor(finalTransfArray[index][0]*trackPts[ti][si].imgPt[li].x + 
                                 finalTransfArray[index][1]*trackPts[ti][si].imgPt[li].y + 
                                 finalTransfArray[index][2]);
		jA = (int) floor(finalTransfArray[index][3]*trackPts[ti][si].imgPt[li].x + 
                                 finalTransfArray[index][4]*trackPts[ti][si].imgPt[li].y + 
                                 finalTransfArray[index][5]);
                
                // check (iA,jA) are inside the image!
		if ( ( iA >= 0) && ( iA < row_max) && ( jA >= 0) && ( jA < col_max)){ 
		  if (is_valid(interpDRG(jA, iA)) ){
                  // calculate ii & jj relative to the image center
		    ii = iA - i_C;
		    jj = jA - j_C;
		    I_e_val = interpDRG(jA,iA) - scaleFactor*trackPts[ti][si].reflectance;
		  }
                  else{
                    I_e_val= 0;
                  }

		  errorArray[index] += abs(I_e_val);
		  
                
		  /*
                  float robustWeight;
                  float b_sqr = 0.0001;
                  if (I_e_val == 0){
		     float tmp = 0.001;
                     robustWeight = sqrt(b_sqr*log(1+(tmp*tmp)/b_sqr))/tmp;
                  }
		  else{
                     robustWeight = sqrt(b_sqr*log(1+(I_e_val*I_e_val)/b_sqr))/fabs(I_e_val);
                  }
                                
                  // We combine the error value with the derivative and
                  // add this to the update equation.
			      	      
		  float weight = spatialWeights(ii+kern_half_width, jj+kern_half_height)*robustWeight;
            
		  */


		  //calculate numerical dirivatives (ii,jj).
                
		  I_x_val = x_deriv(jA,iA)*weight; 
		  I_y_val = y_deriv(jA,iA)*weight; 

		  I_x_I_y = I_x_val*I_y_val;        
		  I_x_sqr = I_x_val*I_x_val;
		  I_y_sqr = I_y_val*I_y_val;

		  // Left hand side
		  lhs(0) += ii * I_x_val * I_e_val;
		  lhs(1) += jj * I_x_val * I_e_val;
		  lhs(2) +=      I_x_val * I_e_val;
		  lhs(3) += ii * I_y_val * I_e_val;
		  lhs(4) += jj * I_y_val * I_e_val;
		  lhs(5) +=      I_y_val * I_e_val;

		  // Right Hand Side UL
		  rhs(0,0) += ii*ii * I_x_sqr;
		  rhs(0,1) += ii*jj * I_x_sqr;
		  rhs(0,2) += ii    * I_x_sqr;
		  rhs(1,1) += jj*jj * I_x_sqr;
		  rhs(1,2) += jj    * I_x_sqr;
		  rhs(2,2) +=         I_x_sqr;

		  // Right Hand Side UR
		  rhs(0,3) += ii*ii * I_x_I_y;
		  rhs(0,4) += ii*jj * I_x_I_y;
		  rhs(0,5) += ii    * I_x_I_y;
		  rhs(1,4) += jj*jj * I_x_I_y;
		  rhs(1,5) += jj    * I_x_I_y;
		  rhs(2,5) +=         I_x_I_y;

		  // Right Hand Side LR
		  rhs(3,3) += ii*ii * I_y_sqr;
		  rhs(3,4) += ii*jj * I_y_sqr;
		  rhs(3,5) += ii    * I_y_sqr;
		  rhs(4,4) += jj*jj * I_y_sqr;
		  rhs(4,5) += jj    * I_y_sqr;
		  rhs(5,5) +=         I_y_sqr;

		}
	      }           
            }// end of if statement: inside image
          }// end of if statement: valid reflectance  
        }// end of for loop over i
      }//end of loop over k

      // Fill in symmetric entries
      rhs(1,0) = rhs(0,1);
      rhs(2,0) = rhs(0,2);
      rhs(2,1) = rhs(1,2);
      rhs(1,3) = rhs(0,4);
      rhs(2,3) = rhs(0,5);
      rhs(2,4) = rhs(1,5);
      rhs(3,0) = rhs(0,3);
      rhs(3,1) = rhs(1,3);
      rhs(3,2) = rhs(2,3);
      rhs(4,0) = rhs(0,4);
      rhs(4,1) = rhs(1,4);
      rhs(4,2) = rhs(2,4);
      rhs(4,3) = rhs(3,4);
      rhs(5,0) = rhs(0,5);
      rhs(5,1) = rhs(1,5);
      rhs(5,2) = rhs(2,5);
      rhs(5,3) = rhs(3,5);
      rhs(5,4) = rhs(4,5);

   
      printRHS(rhs, index, iter);
     
      try {
        solve_symmetric_nocopy(rhs,lhs);
      } catch (ArgumentErr &/*e*/) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
      }


      printIterationValues(errorArray[index], index, iter);

      finalTransfArray[index] += lhs;
  
      iter ++;
 
    }
   
  }//index loop ends here

} 

#if 0
//this is the MP-multi processor version of the above function
void UpdateMatchingParamsLIMA_MP(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
				 ModelParams modelParams, CoregistrationParams coregistrationParams,  
				 vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
				 vector<float> &errorArray )
{

  int numMaxIter = coregistrationParams.maxNumIter;

  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing the scale factor...";
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);
  cout<<"done."<<endl;

  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DRGFilename);

  std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_y_deriv.tif";

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  

  int row_max, col_max;

  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
   
  cout << "image size: rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C = row_max/2;
  int j_C = col_max/2;


  //copy initTransfArray into finalTransfArray.
  for (int index = 0; index < initTransfArray.size(); index++){
    for (int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
  }

   vector<int> ti;
   vector<int> si;
   vector<int> li;
   vector<int> iter;
   vector<int> ii;
   vector<int> jj;
   vector<int> iA;
   vector<int> jA;
   vector<float> I_e_val;
   vector<float> I_x_val;
   vector<float> I_y_val;
   vector<float> I_x_sqr;
   vector<float> I_x_I_y;
   vector<float> I_y_sqr;
   vector<Matrix<float, 6, 6> > rhs;
   vector<Vector<float,6> > lhs; 

   ti.resize(initTransfArray.size());
   si.resize(initTransfArray.size());
   li.resize(initTransfArray.size());
   iter.resize(initTransfArray.size());
   ii.resize(initTransfArray.size());
   jj.resize(initTransfArray.size());
   iA.resize(initTransfArray.size());
   jA.resize(initTransfArray.size());
   I_e_val.resize(initTransfArray.size());
   I_x_val.resize(initTransfArray.size());
   I_y_val.resize(initTransfArray.size());
   I_x_sqr.resize(initTransfArray.size());
   I_x_I_y.resize(initTransfArray.size());
   I_y_sqr.resize(initTransfArray.size());
   rhs.resize(initTransfArray.size());
   lhs.resize(initTransfArray.size());


   int nthreads, tid, chunk; 
   chunk = CHUNKSIZE;
   for (int index = 0; index< initTransfArray.size(); index++){
        iA[index] = 0;
        jA[index] = 0;
        iter[index] = 0;
   }

#pragma omp parallel shared(ti,si, li, iter, ii, jj, iA, jA, I_e_val, I_x_val, I_y_val, I_x_sqr, I_x_I_y, I_y_sqr, rhs, lhs, nthreads,chunk) private(index,tid)
{
  
  tid = omp_get_thread_num();
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
  }
  printf("Thread %d starting...\n",tid);


#pragma omp for schedule(dynamic,chunk)

  for (int index = 0; index < initTransfArray.size(); index++){

    cout << "index = "<< index << endl;

    while( iter[index] <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs[index](i_RHS,j_RHS) = 0.0;
        }
        lhs[index](i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;

      for (ti[index] = 0; ti[index] < trackPts.size(); ti[index]++){

        for (si[index] = 0; si[index] < trackPts[ti[index]].size(); si[index]++){

          if ((trackPts[ti[index]][si[index]].valid ==1) && (trackPts[ti[index]][si[index]].reflectance !=0)){
           
            float weight = trackPts[ti[index]][si[index]].weightRefl;
        
            for (li[index] = 0; li[index] < trackPts[ti[index]][si[index]].LOLAPt.size(); li[index]++){//for each point of a LOLA shot

              if (trackPts[ti[index]][si[index]].LOLAPt[li[index]].s == 3){//center point of a valid shot
    
		iA[index] = (int) floor(finalTransfArray[index][0]*trackPts[ti[index]][si[index]].imgPt[li[index]].x 
                          + finalTransfArray[index][1]*trackPts[ti[index]][si[index]].imgPt[li[index]].y + finalTransfArray[index][2]);
		jA[index] = (int) floor(finalTransfArray[index][3]*trackPts[ti[index]][si[index]].imgPt[li[index]].x 
                          + finalTransfArray[index][4]*trackPts[ti[index]][si[index]].imgPt[li[index]].y + finalTransfArray[index][5]);
                

                // check (iA,jA) are inside the derivative image!
		if ( ( iA[index] >= 0) && ( iA[index] < row_max) && ( jA[index] >= 0) && ( jA[index] < col_max)){

                    	  
                  // calculate ii & jj relative to the image center
		  ii[index] = iA[index] - i_C;
		  jj[index] = jA[index] - j_C;

		  I_e_val[index] = interpDRG(jA[index],iA[index]) - scaleFactor*trackPts[ti[index]][si[index]].reflectance;
		  errorArray[index] += abs(I_e_val[index]);

		  //calculate numerical derivatives (ii,jj).
                  //printf("weight = %f\n", weight);

		  I_x_val[index] = x_deriv(jA[index],iA[index])*weight; 
		  I_y_val[index] = y_deriv(jA[index],iA[index])*weight; 

		  I_x_I_y[index] = I_x_val[index]*I_y_val[index];        
		  I_x_sqr[index] = I_x_val[index]*I_x_val[index];
		  I_y_sqr[index] = I_y_val[index]*I_y_val[index];
                  
                  //printf("I_x_val=%f, I_y_val=%f, I_e_val=%f\n", I_x_val[index],  I_y_val[index], I_e_val[index]);     
		  // Left hand side
		  lhs[index](0) += ii[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](1) += jj[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](2) +=                              I_x_val[index] * I_e_val[index];
		  lhs[index](3) += ii[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](4) += jj[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](5) +=                              I_y_val[index] * I_e_val[index];

		  // Right Hand Side UL
		  rhs[index](0,0) += ii[index]*ii[index] * I_x_sqr[index];
		  rhs[index](0,1) += ii[index]*jj[index] * I_x_sqr[index];
		  rhs[index](0,2) += ii[index]           * I_x_sqr[index];
		  rhs[index](1,1) += jj[index]*jj[index] * I_x_sqr[index];
		  rhs[index](1,2) += jj[index]           * I_x_sqr[index];
		  rhs[index](2,2) +=                       I_x_sqr[index];

		  // Right Hand Side UR
		  rhs[index](0,3) += ii[index]*ii[index] * I_x_I_y[index];
		  rhs[index](0,4) += ii[index]*jj[index] * I_x_I_y[index];
		  rhs[index](0,5) += ii[index]           * I_x_I_y[index];
		  rhs[index](1,4) += jj[index]*jj[index] * I_x_I_y[index];
		  rhs[index](1,5) += jj[index]           * I_x_I_y[index];
		  rhs[index](2,5) +=                       I_x_I_y[index];

		  // Right Hand Side LR
		  rhs[index](3,3) += ii[index]*ii[index] * I_y_sqr[index];
		  rhs[index](3,4) += ii[index]*jj[index] * I_y_sqr[index];
		  rhs[index](3,5) += ii[index]           * I_y_sqr[index];
		  rhs[index](4,4) += jj[index]*jj[index] * I_y_sqr[index];
		  rhs[index](4,5) += jj[index]           * I_y_sqr[index];
		  rhs[index](5,5) +=                       I_y_sqr[index];

		}
	      }           
            }// end of if statement: inside image
          }// end of if statement: valid reflectance  
        }// end of for loop over i
      }//end of loop over k

      // Fill in symmetric entries
      rhs[index](1,0) = rhs[index](0,1);
      rhs[index](2,0) = rhs[index](0,2);
      rhs[index](2,1) = rhs[index](1,2);
      rhs[index](1,3) = rhs[index](0,4);
      rhs[index](2,3) = rhs[index](0,5);
      rhs[index](2,4) = rhs[index](1,5);
      rhs[index](3,0) = rhs[index](0,3);
      rhs[index](3,1) = rhs[index](1,3);
      rhs[index](3,2) = rhs[index](2,3);
      rhs[index](4,0) = rhs[index](0,4);
      rhs[index](4,1) = rhs[index](1,4);
      rhs[index](4,2) = rhs[index](2,4);
      rhs[index](4,3) = rhs[index](3,4);
      rhs[index](5,0) = rhs[index](0,5);
      rhs[index](5,1) = rhs[index](1,5);
      rhs[index](5,2) = rhs[index](2,5);
      rhs[index](5,3) = rhs[index](3,5);
      rhs[index](5,4) = rhs[index](4,5);

      //printf ("before\n");
      //printRHS(rhs[index], index, iter[index]);
      //printLHS_Error(lhs[index], errorArray[index], index, iter[index]);

      try {
        solve_symmetric_nocopy(rhs[index],lhs[index]);
      } catch (ArgumentErr &/*e*/) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
      }


      printLHS_Error(lhs[index], errorArray[index], index, iter[index]);

      finalTransfArray[index] += lhs[index];
  
      iter[index] ++;
 
    }
   
  }//index loop ends here

 }//openMP


} 
#endif


//this is the MP-multi processor version of the above function
void UpdateMatchingParamsLIDEM_MP(vector<vector<LOLAShot> > &trackPts, string DEMFilename,  
				  ModelParams modelParams,  CoregistrationParams coregistrationParams, 
				  vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
				  vector<float> &errorArray )
{
#if 0
  int numMaxIter = coregistrationParams.maxNumIter;

  DiskImageView<PixelGray<float> >   DEM(DEMFilename);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, DEMFilename);

  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<float> >   interpDEM = interpolate(edge_extend(DEM.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DEMFilename);

  std::string xDerivFilename = "../results" + prefix_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_from_filename(temp) + "_y_deriv.tif";

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DEM),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DEM),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  

  int row_max, col_max;

  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
   
  cout << "image size: rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C = row_max/2;
  int j_C = col_max/2;


  //copy initTransfArray into finalTransfArray.
  for (int index = 0; index < initTransfArray.size(); index++){
    for (int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
  }

  
   vector<int> ti;
   vector<int> si;
   vector<int> li;
   vector<int> iter;
   vector<int> ii;
   vector<int> jj;
   vector<int> iA;
   vector<int> jA;
   vector<float> I_e_val;
   vector<float> I_x_val;
   vector<float> I_y_val;
   vector<float> I_x_sqr;
   vector<float> I_x_I_y;
   vector<float> I_y_sqr;
   vector<Matrix<float, 6, 6> > rhs;
   vector<Vector<float,6> > lhs; 

   ti.resize(initTransfArray.size());
   si.resize(initTransfArray.size());
   li.resize(initTransfArray.size());
   iter.resize(initTransfArray.size());
   ii.resize(initTransfArray.size());
   jj.resize(initTransfArray.size());
   iA.resize(initTransfArray.size());
   jA.resize(initTransfArray.size());
   I_e_val.resize(initTransfArray.size());
   I_x_val.resize(initTransfArray.size());
   I_y_val.resize(initTransfArray.size());
   I_x_sqr.resize(initTransfArray.size());
   I_x_I_y.resize(initTransfArray.size());
   I_y_sqr.resize(initTransfArray.size());
   rhs.resize(initTransfArray.size());
   lhs.resize(initTransfArray.size());

   int index;
   int nthreads, tid, chunk;
     
   chunk = CHUNKSIZE;
   for (int index = 0; index< initTransfArray.size(); index++){
        iA[index] = 0;
        jA[index] = 0;
        iter[index] = 0;
   }

#pragma omp parallel shared(ti, si, li, iter, ii, jj, iA, jA, I_e_val, I_x_val, I_y_val, I_x_sqr, I_x_I_y, I_y_sqr, rhs, lhs, nthreads,chunk) private(index,tid)
{
  
  tid = omp_get_thread_num();
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
  }
  printf("Thread %d starting...\n",tid);


#pragma omp for schedule(dynamic,chunk)

  for (int index = 0; index < initTransfArray.size(); index++){
    //for (index = 124; index < 128; index++){
 
    while( iter[index] <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs[index](i_RHS,j_RHS) = 0.0;
        }
        lhs[index](i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;

      for (ti[index] = 0; ti[index] < trackPts.size(); ti[index]++){

        for (si[index] = 0; si[index] < trackPts[ti[index]].size(); si[index]++){

          if (trackPts[ti[index]][si[index]].valid ==1){
           
            //weight = trackPts[k][i].weight_prd;
        
            for (li[index] = 0; li[index] < trackPts[ti[index]][si[index]].LOLAPt.size(); li[index]++){//for each point of a LOLA shot

     
		iA[index] = (int) floor(finalTransfArray[index][0]*trackPts[ti[index]][si[index]].DEMPt[li[index]].x 
                          + finalTransfArray[index][1]*trackPts[ti[index]][si[index]].DEMPt[li[index]].y + finalTransfArray[index][2]);
		jA[index] = (int) floor(finalTransfArray[index][3]*trackPts[ti[index]][si[index]].DEMPt[li[index]].x 
                          + finalTransfArray[index][4]*trackPts[ti[index]][si[index]].DEMPt[li[index]].y + finalTransfArray[index][5]);
                
                // check (iA,jA) are inside the image!
		if ( ( iA[index] >= 0) && ( iA[index] < row_max) && ( jA[index] >= 0) && ( jA[index] < col_max)){
        	  
                  // calculate ii & jj relative to the image center
		  ii[index] = iA[index] - i_C;
		  jj[index] = jA[index] - j_C;
                 
		  I_e_val[index] = interpDEM(jA[index],iA[index]) - trackPts[ti[index]][si[index]].LOLAPt[li[index]].coords[2];
               
		  errorArray[index] += abs(I_e_val[index]);

		  //calculate numerical dirivatives (ii,jj).
		  I_x_val[index] = x_deriv(jA[index],iA[index]); 
		  I_y_val[index] = y_deriv(jA[index],iA[index]); 

		  I_x_I_y[index] = I_x_val[index]*I_y_val[index];        
		  I_x_sqr[index] = I_x_val[index]*I_x_val[index];
		  I_y_sqr[index] = I_y_val[index]*I_y_val[index];

		  // Left hand side
		  lhs[index](0) += ii[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](1) += jj[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](2) +=                              I_x_val[index] * I_e_val[index];
		  lhs[index](3) += ii[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](4) += jj[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](5) +=                              I_y_val[index] * I_e_val[index];

		  // Right Hand Side UL
		  rhs[index](0,0) += ii[index]*ii[index] * I_x_sqr[index];
		  rhs[index](0,1) += ii[index]*jj[index] * I_x_sqr[index];
		  rhs[index](0,2) += ii[index]           * I_x_sqr[index];
		  rhs[index](1,1) += jj[index]*jj[index] * I_x_sqr[index];
		  rhs[index](1,2) += jj[index]           * I_x_sqr[index];
		  rhs[index](2,2) +=                       I_x_sqr[index];

		  // Right Hand Side UR
		  rhs[index](0,3) += ii[index]*ii[index] * I_x_I_y[index];
		  rhs[index](0,4) += ii[index]*jj[index] * I_x_I_y[index];
		  rhs[index](0,5) += ii[index]           * I_x_I_y[index];
		  rhs[index](1,4) += jj[index]*jj[index] * I_x_I_y[index];
		  rhs[index](1,5) += jj[index]           * I_x_I_y[index];
		  rhs[index](2,5) +=                       I_x_I_y[index];

		  // Right Hand Side LR
		  rhs[index](3,3) += ii[index]*ii[index] * I_y_sqr[index];
		  rhs[index](3,4) += ii[index]*jj[index] * I_y_sqr[index];
		  rhs[index](3,5) += ii[index]           * I_y_sqr[index];
		  rhs[index](4,4) += jj[index]*jj[index] * I_y_sqr[index];
		  rhs[index](4,5) += jj[index]           * I_y_sqr[index];
		  rhs[index](5,5) +=                       I_y_sqr[index];

		}//end of if statemenet: check (iA,jA) are inside the image! 
            }// end of for loop over li
          }// end of if statement: valid track and reflectance  
        }// end of for loop over si
      }//end of loop over ti

      // Fill in symmetric entries
      rhs[index](1,0) = rhs[index](0,1);
      rhs[index](2,0) = rhs[index](0,2);
      rhs[index](2,1) = rhs[index](1,2);
      rhs[index](1,3) = rhs[index](0,4);
      rhs[index](2,3) = rhs[index](0,5);
      rhs[index](2,4) = rhs[index](1,5);
      rhs[index](3,0) = rhs[index](0,3);
      rhs[index](3,1) = rhs[index](1,3);
      rhs[index](3,2) = rhs[index](2,3);
      rhs[index](4,0) = rhs[index](0,4);
      rhs[index](4,1) = rhs[index](1,4);
      rhs[index](4,2) = rhs[index](2,4);
      rhs[index](4,3) = rhs[index](3,4);
      rhs[index](5,0) = rhs[index](0,5);
      rhs[index](5,1) = rhs[index](1,5);
      rhs[index](5,2) = rhs[index](2,5);
      rhs[index](5,3) = rhs[index](3,5);
      rhs[index](5,4) = rhs[index](4,5);

   
      //printRHS(rhs[index], index, iter[index]);
     
      try {
        solve_symmetric_nocopy(rhs[index],lhs[index]);
      } catch (ArgumentErr &/*e*/) {
        //             std::cout << "Error @ " << x << " " << y << "\n";
        //             std::cout << "Exception caught: " << e.what() << "\n";
        //             std::cout << "PRERHS: " << pre_rhs << "\n";
        //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
        //             std::cout << "RHS: " << rhs << "\n";
        //             std::cout << "LHS: " << lhs << "\n\n";
        //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
        //             exit(0);
      }


      //printLHS_Error(lhs[index], errorArray[index], index, iter[index]);

      finalTransfArray[index] += lhs[index];
  
      iter[index] ++;
 
    }
   
  }//index loop ends here

 }//openMP
#endif

} 
#endif

