// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <valarray>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

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

//computes the translation between the foreground and background pixels
Vector3 ComputeDEMTranslation( const vector<Vector3>& features, 
                               const vector<Vector3>& reference ) {
  Vector3 translation(0,0,0);
  
  unsigned int numValidMatches = 0;
  for( unsigned int i = 0; i < features.size(); ++i ){
    if( norm_1(reference[i]) > 0 ){
      translation += reference[i] - features[i];
      ++numValidMatches;
    }
  }
  return translation/numValidMatches;
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

//compute the matching error vector in the form of a standard deviation 
valarray<float> ComputeMatchingError( const vector<Vector3>& model, 
                                      const vector<Vector3>& reference ) {
  if( model.size() != reference.size() ){
    vw_throw( ArgumentErr() << 
    "The vectors passed to ComputeMatchingError() have different sizes." ); 
  }

  valarray<float> errors( model.size() );
  for( unsigned int i = 0; i < model.size(); ++i ){
    errors[i] = norm_2( model[i] - reference[i] );
  }
  return errors;
}

// This used to return a Vector3 with errors in each dimension, now just returns total error.
//compute the matching error vector in the form of individual errors
float ComputeMatchingError3D( const vector<Vector3>& model, 
                              const vector<Vector3>& reference ) {
  valarray<float> errors = ComputeMatchingError( model, reference );
  return errors.sum()/errors.size();
}

Matrix<float, 3, 3> ComputeDEMRotation( const vector<Vector3>& features, 
                                        const vector<Vector3>& matches,
                                        const Vector3&         featureCenter,
                                        const Vector3&         matchCenter ) {
  Matrix<float, 3, 3> rotation;
 
  Matrix<float,3,3> A;
  Matrix<float,3,3> U;
  Vector<float,3>   s;
  Matrix<float,3,3> V;
  
  A.set_zero();
  
  vw_out(vw::InfoMessage, "icp") << "  F_center: " << featureCenter << endl;
  vw_out(vw::InfoMessage, "icp") << "  M_center: " << matchCenter   << endl;
  
  for( unsigned int i = 0; i < features.size(); ++i ){
    if( norm_1(matches[i]) > 0 ){ //ignore invalid matches
      A += outer_prod( matches[i] - matchCenter, features[i] - featureCenter );
    }
  }
  
  svd(A, U, s, V);

  // use Kabsch Algorithm to form the rotation matrix
  if( det(A) < 0 ){
    Matrix3x3 sign_id = identity_matrix(3);
    sign_id(2,2) = -1;
    rotation = U*sign_id*V;
  } 
  else { rotation = U*V; }
 
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
//applies a 3D rotation and transform to a DEM
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 featureCenter, Vector3 translation, Matrix<float,3,3> rotation)
{

  for (unsigned int i=0; i < featureArray.size(); i++){
    featureArray[i] = rotation*(featureArray[i]-featureCenter) + featureCenter+translation; 
    //featureArray[i] = rotation*(featureArray[i]) + translation; 
  }

 
}






