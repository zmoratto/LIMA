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

//computes the translation between the foreground and background pixels
void ComputeDEMTranslation(vector<Vector3> featureArray, vector<Vector3> matchArray, Vector3 &translation )
{

  //Vector3 translation;
  translation[0] = 0;
  translation[1] = 0;
  translation[2] = 0;

  int count = 0;
  
  for (int i = 0; i < featureArray.size(); i++){
       Vector3 fore_lonlat3 = featureArray[i];
       Vector3 back_lonlat3 = matchArray[i];
       translation[0] = translation[0] + fore_lonlat3(0) - back_lonlat3(0);
       translation[1] = translation[1] + fore_lonlat3(1) - back_lonlat3(1);
       translation[2] = translation[2] + fore_lonlat3(2) - back_lonlat3(2);
       count++;
  }

  translation[0] = translation[0]/count;
  translation[1] = translation[1]/count;
  translation[2] = translation[2]/count;

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
float ComputeMatchingError(vector<Vector3> featureArray, vector<Vector3> matchArray, vector<float> &errorArray)
{


   float matchError = 0.0;
   for (int i = 0; i < featureArray.size(); i++){
       float overallDist = 0.0;
       float dist;
       for (int j = 0; j < 3; j++){
         dist = featureArray[i](j) - matchArray[i](j);
         //cout<<j<<" "<<featureArray[i](j)<<" "<<matchArray[i](j)<<" "<<dist<<endl;
	 overallDist = overallDist + dist*dist; 
       }
       //cout<<endl;
       errorArray[i]=sqrt(overallDist);
       matchError = matchError + errorArray[i];
   }
   return matchError/featureArray.size();
  
}

void ComputeDEMRotation(vector<Vector3> featureArray, vector<Vector3> matchArray, 
                        Vector3 translation, Matrix<float, 3, 3> &rotation)
{

  
  //Matrix<float, 3, 3> rotation;
  Matrix<float,3,3> A;
  Matrix<float,3,3> U;
  Vector<float,3>   s;
  Matrix<float,3,3> V;

  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      A[0][0] = 0.0;
    }
  }


  for (int i = 0; i < featureArray.size(); i++){     
       
       Vector3 fore_lonlat3 = featureArray[i];
       Vector3 back_lonlat3 = matchArray[i];

       A[0][0] = A[0][0] + back_lonlat3[0]*(fore_lonlat3[0] - translation[0]);
       A[0][1] = A[0][1] + back_lonlat3[1]*(fore_lonlat3[0] - translation[0]);
       A[0][2] = A[0][2] + back_lonlat3[2]*(fore_lonlat3[0] - translation[0]);
       A[1][0] = A[1][0] + back_lonlat3[0]*(fore_lonlat3[1] - translation[1]);
       A[1][1] = A[1][1] + back_lonlat3[1]*(fore_lonlat3[1] - translation[1]);
       A[1][2] = A[1][2] + back_lonlat3[2]*(fore_lonlat3[1] - translation[1]);
       A[2][0] = A[2][0] + back_lonlat3[0]*(fore_lonlat3[2] - translation[2]);
       A[2][1] = A[2][1] + back_lonlat3[1]*(fore_lonlat3[2] - translation[2]);
       A[2][2] = A[2][2] + back_lonlat3[2]*(fore_lonlat3[2] - translation[2]);
  }
  

  //printf("correlation: A\n");
  //PrintMatrix(A);
  
  //extract the rotation matrix
  //void svd( AMatrixT const& A, UMatrixT &U, SingularValuesT &s, VTMatrixT &VT );
  svd(A, U, s, V);
  

  //printf("A = UWV^T\n");
  //printf("W\n");
  //PrintMatrix(A);
 
  //printf("U\n");
  //PrintMatrix(U);
 
  //printf("V\n");
  //PrintMatrix(V);

  //printf("SVD done\n");

  Matrix<float,3,3> VT = transpose(V);
  //PrintMatrix(VT);

  rotation = U*VT;

  //printf("R\n");
  //PrintMatrix(rotation);

  //return rotation;

}

//applies a 3D rotation and transform to a DEM
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 translation, Matrix<float,3,3> rotation)
{

  for (int i=0; i < featureArray.size(); i++){
    featureArray[i] = rotation*featureArray[i] + translation; 
  }

 
}



