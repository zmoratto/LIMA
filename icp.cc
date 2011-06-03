// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <valarray>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

#include "icp.h"
#include "coregister.h"

Vector2 fore_2_back_lonlat(Vector2 fore_lon_lat)
{ 
  //change into USGS coords
  float usgs_2_lonlat = 180/(3.14159265*3396190);
  Vector2 back_lon_lat;
  back_lon_lat(0) = fore_lon_lat(0)-180;  
  back_lon_lat(1) = fore_lon_lat(1);  
  back_lon_lat = back_lon_lat/usgs_2_lonlat;
  return back_lon_lat;
} 

Vector2 back_2_fore_lonlat(Vector2 back_lon_lat)
{
  float usgs_2_lonlat = 180/(3.14159265*3396190);
  Vector2 fore_lon_lat;
  fore_lon_lat(0) = back_lon_lat(0)*usgs_2_lonlat+180;
  fore_lon_lat(1) = back_lon_lat(1)*usgs_2_lonlat; 
  return fore_lon_lat;
}

//computes the translation between the foreground and background pixels
Vector3 ComputeDEMTranslation(const vector<Vector3>& featureArray, 
	                      const vector<Vector3>& matchArray)
{
  Vector3 translation(0,0,0);
  
  int numValidMatches = 0;
  for (unsigned int i = 0; i < featureArray.size(); i++){
    if ((matchArray[i][0]!=0) && (matchArray[i][1]!=0) && (matchArray[i][2]!=0)){
       translation += matchArray[i] - featureArray[i];
       numValidMatches++;
    }
  }

  return translation / numValidMatches;//featureArray.size();
}


void PrintMatrix(Matrix<float, 3, 3> &A)
{
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      printf("%f ", A[i][j]);
    }
    printf("\n");
  }
}

//compute the matching error vector and overall average error 
valarray<float> 
ComputeMatchingError(
	const vector<Vector3>& featureArray, 
	const vector<Vector3>& matchArray)
{
valarray<float> errorArray( featureArray.size() );
for (unsigned int i = 0; i < featureArray.size(); i++)
	{
	float overallDist(0.0);
	float dist(0.0);
	for (int j = 0; j < 3; j++)
		{
		dist = featureArray[i](j) - matchArray[i](j);
		//cout<<j<<" "<<featureArray[i](j)<<" "<<matchArray[i](j)<<" "<<dist<<endl;
		overallDist += dist*dist; 
		}
	errorArray[i]=sqrt(overallDist);
	//cout<<"error= "<<errorArray[i]<<endl;
	}
return errorArray;
}

Matrix<float, 3, 3> ComputeDEMRotation(const vector<Vector3>& featureArray, 
	                               const vector<Vector3>& matchArray)
	// Vector3 translation, <-- This used to be passed in, but not used???
{

  Matrix<float, 3, 3> rotation;
  
  Matrix<float,3,3> A;
  Matrix<float,3,3> U;
  Vector<float,3>   s;
  Matrix<float,3,3> V;
  
  A.set_zero(); // Fill with zeros, maybe what the below meant?

  //compute the centroids
  Vector3 featureCenter;
  Vector3 matchCenter;
  int numValidMatches = 0;
  for (unsigned int i = 0; i < featureArray.size(); i++){ 
    if ((matchArray[i][0]!=0) && (matchArray[i][1]!=0) && (matchArray[i][2]!=0)){//ignore invalid matches
       featureCenter += featureArray[i];
       matchCenter += matchArray[i];
       numValidMatches++;
    }
  }
  featureCenter /= numValidMatches;
  matchCenter   /= numValidMatches;
  
  vw_out(vw::InfoMessage, "icp") << "  F_center: " << featureCenter << endl;
  vw_out(vw::InfoMessage, "icp") << "  M_center: " << matchCenter   << endl;
  
  for (unsigned int i = 0; i < featureArray.size(); i++){     
    Vector3 feature = featureArray[i];
    Vector3 match = matchArray[i];
    
    if ((matchArray[i][0]!=0) && (matchArray[i][1]!=0) && (matchArray[i][2]!=0)){//ignore invalid matches
    
      A[0][0] = A[0][0] + (match[0]-matchCenter[0])*(feature[0] - featureCenter[0]);
      A[1][0] = A[1][0] + (match[1]-matchCenter[1])*(feature[0] - featureCenter[0]);
      A[2][0] = A[2][0] + (match[2]-matchCenter[2])*(feature[0] - featureCenter[0]);
      
      A[0][1] = A[0][1] + (match[0]-matchCenter[0])*(feature[1] - featureCenter[1]);
      A[1][1] = A[1][1] + (match[1]-matchCenter[1])*(feature[1] - featureCenter[1]);
      A[2][1] = A[2][1] + (match[2]-matchCenter[2])*(feature[1] - featureCenter[1]);
      
      A[0][2] = A[0][2] + (match[0]-matchCenter[0])*(feature[2] - featureCenter[2]);
      A[1][2] = A[1][2] + (match[1]-matchCenter[1])*(feature[2] - featureCenter[2]);
      A[2][2] = A[2][2] + (match[2]-matchCenter[2])*(feature[2] - featureCenter[2]);
    }
  }
  
svd(A, U, s, V);
 
/* Never used? 
Matrix<float,3,3> VT = transpose(V);
*/

rotation = U*V;

/* Never used?
Matrix<float,3,3> id = rotation*transpose(rotation); 
*/

return rotation;
}

//applies a 3D rotation and transform to a DEM
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 translation, Matrix<float,3,3> rotation)
{


  Vector3 featureCenter;
  for (unsigned int i = 0; i < featureArray.size(); i++){  
      featureCenter = featureCenter + featureArray[i];
  }
  featureCenter = featureCenter/featureArray.size();

  for (unsigned int i=0; i < featureArray.size(); i++){
    featureArray[i] = rotation*(featureArray[i]-featureCenter) + featureCenter+translation; 
  }

 
}






