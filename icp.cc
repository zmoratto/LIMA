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

#include "util.h"
#include "icp.h"
#include "coregister.h"

//computes the translation T = avg(reference-match)
Vector3 ComputeDEMTranslation( const vector<Vector3>& match, 
                               const vector<Vector3>& reference ) {
  Vector3 translation(0,0,0);
  
  unsigned int numValidMatches = 0;
  for( unsigned int i = 0; i < match.size(); ++i ){
    if( norm_1(reference[i]) > 0 ){
      translation += reference[i] - match[i];
      ++numValidMatches;
    }
  }
  return translation/numValidMatches;
}


//computes the translation T = avg(reference-rotation*match)
Vector3 ComputeDEMTranslation( const vector<Vector3>& match, 
                               const vector<Vector3>& reference,
                               const Matrix<float, 3, 3>& rotation) {
  Vector3 translation(0,0,0);
  
  unsigned int numValidMatches = 0;
  
  //PrintMatrix(rotation);
  
  for( unsigned int i = 0; i < match.size(); ++i ){
    if( norm_1(reference[i]) > 0 ){
      translation += reference[i] - rotation*match[i];
      ++numValidMatches;
    }
  }
  return translation/numValidMatches;
}

void PrintMatrix(const Matrix<float, 3, 3> &A)
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
//compute a rotation matrix between two point clouds
Matrix<float, 3, 3> ComputeDEMRotation( const vector<Vector3>& match, 
                                        const vector<Vector3>& reference,
                                        const Vector3&         matchCenter,
                                        const Vector3&         referenceCenter ) {
  Matrix<float, 3, 3> rotation;
 
  Matrix<float,3,3> A;
  Matrix<float,3,3> U;
  Vector<float,3>   s;
  Matrix<float,3,3> V;
  
  A.set_zero();
  
  vw_out(vw::InfoMessage, "icp") << "  Match_center: " << matchCenter << endl;
  vw_out(vw::InfoMessage, "icp") << "  Reference_center: " << referenceCenter   << endl;
  
  for( unsigned int i = 0; i < match.size(); ++i ){
    if( norm_1(reference[i]) > 0 ){ //ignore invalid matches
   
      A += outer_prod( reference[i] - referenceCenter, match[i] - matchCenter );
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

//computes the rotation matrix R
//such that reference = R*match
Matrix<float, 3, 3> ComputeDEMRotation( const vector<Vector3>& match, 
                                        const vector<Vector3>& reference,
                                        const Vector3&         translation) {
  Matrix<float, 3, 3> rotation;
 
  Matrix<float,3,3> A;
  Matrix<float,3,3> U;
  Vector<float,3>   s;
  Matrix<float,3,3> V;
  
  A.set_zero();
  
  vw_out(vw::InfoMessage, "icp") << "  translation: " << translation << endl;

  
  for( unsigned int i = 0; i < match.size(); ++i ){
    if( norm_1(reference[i]) > 0 ){ //ignore invalid matches
      A += outer_prod( reference[i], match[i] + translation);
      //correction to the order in which the matrix is filled
      //A += outer_prod( match[i] + translation, reference[i]);
    }
  }
  
  svd(A, U, s, V);

  // use Kabsch Algorithm to form the rotation matrix
  if( det(A) < 0 ){
    Matrix3x3 sign_id = identity_matrix(3);
    sign_id(2,2) = -1;
    rotation = U*sign_id*V;
    //cout<<"DET NEGATIVE"<<endl;
  } 
  else { rotation = U*V; 
    //cout<<"DET POSITIVE"<<endl;
  }
 
  return rotation;
}

//applies a 3D rotation and translation to a set of 3D points
//extracts the centroid of the original set of 3D points
void TransformFeatures( vector<Vector3>&         featureArray, 
                        const Vector3&           translation, 
                        const Matrix<float,3,3>& rotation ){
  const Vector3 featureCenter = find_centroid( featureArray );
  TransformFeatures( featureArray, featureCenter, translation, rotation );
}

//applies a 3D rotation and translation to a set of 3D points
//assumes a precomputed centroid of the 3D point set
void TransformFeatures(       vector<Vector3>&   featureArray, 
                        const Vector3&           featureCenter, 
                        const Vector3&           translation, 
                        const Matrix<float,3,3>& rotation) {
  for( unsigned int i=0; i < featureArray.size(); ++i ){
    featureArray[i] = rotation*(featureArray[i]-featureCenter) + featureCenter + translation;      
  }
}

//used to check the quality of matches for debug purposes
void ICP_Report(std::vector<vw::Vector3>       origFeatureXYZArray, 
		std::vector<vw::Vector3>       featureXYZArray, 
		std::vector<vw::Vector3>       referenceXYZArray, 
		const vw::cartography::GeoReference& backDEMGeo, 
		const vw::cartography::GeoReference& foreDEMGeo,
		vw::Vector3&                   center, 
		vw::Vector3&                   translation, 
		vw::Matrix<float, 3, 3>&       rotation,
                float matchingError)
{
 
    //print the feature and reference arrays for debug - START 
    for (int index = 0; index<featureXYZArray.size(); index++){
    
      Vector3 transfFeatureXYZ = rotation*(origFeatureXYZArray[index]/*-center*/) /*+ center */+ translation; 

      Vector3 origFeatureLonLatRad = foreDEMGeo.datum().cartesian_to_geodetic(origFeatureXYZArray[index]);
      Vector3 featureLonLatRad = foreDEMGeo.datum().cartesian_to_geodetic(featureXYZArray[index]);
      Vector3 transFeatureLonLatRad = foreDEMGeo.datum().cartesian_to_geodetic(transfFeatureXYZ);
      Vector3 referenceLonLatRad = foreDEMGeo.datum().cartesian_to_geodetic(referenceXYZArray[index]);
      
      cout<<"ICP report: origFeature_xyz: "<<origFeatureXYZArray[index]<<endl;
      cout<<"ICP report: transfFeature_xyz: "<<transfFeatureXYZ<<endl;
      cout<<"ICP report: feature_xyz: "<<featureXYZArray[index]<<endl;
      cout<<"ICP report:reference_xyz: "<<referenceXYZArray[index]<<endl;
      //Vector3 featureLonLatRad = backDEMGeo.datum().cartesian_to_geodetic(featureXYZArray[index]);
      //Vector3 referenceLonLatRad = backDEMGeo.datum().cartesian_to_geodetic(referenceXYZArray[index]);
      cout<<"ICP report:origFeature_llr: "<<origFeatureLonLatRad<<endl;
      cout<<"ICP report:transfFeature_llr: "<<transFeatureLonLatRad<<endl;
      cout<<"ICP report:feature_llr: "<<featureLonLatRad<<endl;
      cout<<"ICP report:reference_llr: "<<referenceLonLatRad<<endl;
      cout<<"ICP report: matching error: "<<matchingError<<endl;
      cout<<"-----------------------------------------------------"<<endl;
    }
    //print the feature and reference arrays for debug - END
}




