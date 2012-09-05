// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifndef ICP_H
#define ICP_H

#include <string>
#include <valarray>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include "coregister.h"
#include "util.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

vw::Vector3 ComputeDEMTranslation( const std::vector<vw::Vector3>& features, 
                                   const std::vector<vw::Vector3>& reference );

vw::Vector3 ComputeDEMTranslation( const vector<Vector3>& match, 
                                   const vector<Vector3>& reference,
                                   const Matrix<float, 3, 3> &rotation,
                                   const Vector3&              matchCenter );

void PrintMatrix(const Matrix<float, 3, 3> &A);

//compute the matching error vector and overall average error 
void ComputeMatchingError(       std::valarray<float>& errors,
                           const std::vector<vw::Vector3>& model, 
                           const std::vector<vw::Vector3>& reference );

float ComputeMatchingError3D( const std::vector<vw::Vector3>& model, 
                              const std::vector<vw::Vector3>& reference );

Matrix<float, 3, 3> ComputeDEMRotation( const std::vector<vw::Vector3>& features, 
                                        const std::vector<vw::Vector3>& matches,
                                        const vw::Vector3&              featureCenter,
                                        const vw::Vector3&              matchCenter );

inline Matrix<float, 3, 3> ComputeDEMRotation( const vector<Vector3>& features, 
                                               const vector<Vector3>& matches ) {
  const Vector3 featureCenter = find_centroid( features );
  const Vector3 matchCenter = find_centroid( matches );
  return ComputeDEMRotation( features, matches, featureCenter, matchCenter);
}

Matrix<float, 3, 3> ComputeDEMRotation( const vector<Vector3>& match, 
                                        const vector<Vector3>& reference,
                                        const Vector3&         translation);

//applies a 3D rotation and transform to a DEM
void TransformFeatures(       std::vector<vw::Vector3>& featureArray, 
                        const vw::Vector3&              translation, 
                        const vw::Matrix<float,3,3>&    rotation );

void TransformFeatures(       std::vector<vw::Vector3>& featureArray, 
                        const vw::Vector3&              featureCenter, 
                        const vw::Vector3&              translation, 
                        const vw::Matrix<float,3,3>&    rotation );

//compute a set of features from the fore image.
template <class ViewT>
std::vector<vw::Vector3> GetFeatures( const vw::ImageViewBase<ViewT>&      im,
                                      const vw::cartography::GeoReference& georef,
                                      const vw::Vector2&                   samplingStep, 
                                      const vw::Vector2&                   delta_lonlat,
                                      const float                          noDataVal ) {
  //cout<<"Get Features..."<<endl;
  std::vector<vw::Vector3> features;
  
  for( int j = 0; j < im.impl().rows(); j += samplingStep.y() ){
    for( int i = 0; i < im.impl().cols(); i += samplingStep.x() ){
      //determine the location of the corresponding background point for each valid fore point
      if( (isnan(im.impl()(i,j))!=FP_NAN) && (im.impl()(i,j) != noDataVal) ){
        vw::Vector2 lonlat = georef.pixel_to_lonlat( vw::Vector2(i,j) );
        lonlat += delta_lonlat;
        
        if (lonlat(0) < -180){
	  lonlat(0)=360+lonlat(0);
        }
	
        if (lonlat(0) > 180){
	  lonlat(0)=360-lonlat(0);
        }
	
	vw::Vector3 lonlatrad( lonlat.x(), lonlat.y(), im.impl()(i,j) );
	vw::Vector3 xyz = georef.datum().geodetic_to_cartesian( lonlatrad );
	
	features.push_back( xyz );
      }
    }
  }
  //cout<<"numFeatures="<<featureArray.size()<<endl;
  return features;
};

template <class ViewT>
std::vector<vw::Vector3> GetFeatures( const vw::ImageViewBase<ViewT>&      foreImg,
                                      const vw::cartography::GeoReference& foreGeo,
                                      const vw::ImageViewBase<ViewT>&      backImg, 
                                      const vw::cartography::GeoReference& backGeo, 
                                      const vw::Vector2&                   samplingStep, 
                                      const vw::Vector2&                   delta_lonlat,
                                      const float                          noDataVal ) {
  // Once upon a time this function may have used the backImg, but it doesn't anymore.
  // This interface is kept around for historic reasons.
  return GetFeatures( foreImg, foreGeo, samplingStep, delta_lonlat, noDataVal );
}

//find the closest background points to the fore points.
//used in ICP_DEM_2_DEM
template <class ViewT>
vector<bool> FindMatches( const std::vector<vw::Vector3>&      features, 
                          const vw::ImageViewBase<ViewT>&      backImg, 
                          const vw::cartography::GeoReference& backGeo, 
                          const vw::cartography::GeoReference& foreGeo, 
                                std::vector<vw::Vector3>&      referenceArray, 
                          const vw::Vector2&                   matchWindowHalfSize, 
                          const float                          noValData) {
  vector<bool> mask( referenceArray.size(), true );
  
  vw::ImageViewRef<typename ViewT::pixel_type> interpBackImg = 
        interpolate(edge_extend(backImg.impl(), vw::ConstantEdgeExtension()), 
                    vw::BilinearInterpolation()); 

  for( unsigned int i = 0; i < features.size(); ++i ){
    //revert fore from cartesian to spherical coordinates
    vw::Vector3 fore_lonlat3 = foreGeo.datum().cartesian_to_geodetic( features[i] );

    vw::Vector2 backPix = backGeo.lonlat_to_pixel( vw::Vector2( fore_lonlat3.x(), 
                                                                fore_lonlat3.y() ) );

    vw::Vector3 back_xyz; //invalid point
    back_xyz[0]=0.0;
    back_xyz[1]=0.0;
    back_xyz[2]=0.0;

    float minDistance = std::numeric_limits<float>::max();
    vw::Vector<int,2> start = backPix - matchWindowHalfSize;
    vw::Vector<int,2> limit = backPix + matchWindowHalfSize + vw::Vector<int,2>( 1, 1 );

    for( int k = start.y(); k < limit.y(); ++k ){
      for( int l = start.x(); l < limit.x(); ++l ){
        float interpBack = interpBackImg(l,k);
        if( interpBack == noValData ){
           continue; 
        }

	//cout<<"fore="<<features[i][0]<<", "<<features[i][1]<<", "<<features[i][2]<<", back="<<interpBack<<endl;
        //vw::Vector2 bPixel = backGeo.pixel_to_lonlat( vw::Vector2(l,k) );
	vw::Vector2 back_lonlat = backGeo.pixel_to_lonlat( vw::Vector2(l,k) );
	//std::cout<<bPixel(0)<<", "<<bPixel(1)<<std::endl;  

        //revert to fore lon lat system
        vw::Vector3 t_back_lonlat3( back_lonlat.x(), back_lonlat.y(), interpBack );
        //cout<<"fore_lonlat3="<<fore_lonlat3<<endl;
        //cout<<"back_lonlat3"<<t_back_lonlat3<<endl;
     
        //transform into xyz coordinates of the foregound image.
        vw::Vector3 t_back_xyz= foreGeo.datum().geodetic_to_cartesian(t_back_lonlat3);
        //cout<<"fore_xyz="<<features[i]<<endl;
        //cout<<"back_xyz="<<t_back_xyz<<endl; 
       
        float distance = norm_2( t_back_xyz - features[i] );
	
        if( distance < minDistance ){
          minDistance = distance;
          back_xyz = t_back_xyz;
        }
      }
    }
    referenceArray[i]=back_xyz;
    if( back_xyz == vw::Vector3(0,0,0) ){
      mask[i] = false;
    }
  }
  return mask;
};

//determines best DEM points that match to LOLA track points
//used in ICP and in error calculations.
//the matching vector is returned in cartesian coordinates.
//the match is done by minimizing the Euclidian distance in the 3D cartesian space.
//TODO: add hor and ver pixel offset together with altitude
template <class ViewT>
//void FindMatchesFromDEM( const std::vector<vw::Vector3>&      xyzModel, 
vector<bool> FindMatchesFromDEM( const std::vector<vw::Vector3>&      xyzModel, 
                                 const std::vector<vw::Vector3>&      llrModel,
                                 const vw::ImageViewBase<ViewT>&      DEM, 
                                 const vw::cartography::GeoReference& DEMGeo, 
                                       std::vector<vw::Vector3>&      xyzMatch, 
                                 const vw::Vector3&                   translation, 
                                 const vw::Matrix<float,3,3>&         rotation, 
                                 const vw::Vector3&                   xyzMatchCenter,
                                 const double                         noDEMVal, 
                                 const vw::Vector2&                   matchWindowHalfSize ) {
  vector<bool> mask( xyzMatch.size(), true );
  
  // cout<<"XYZ_MATCH_CENTER="<<xyzMatchCenter<<endl;
  // cout<<DEMGeo<<endl;  

  vw::BBox2 DEM_bbox = bounding_box( DEM );
  for( unsigned int index = 0; index < xyzModel.size(); ++index ){
    const vw::Vector2 DEM_pix = DEMGeo.lonlat_to_pixel( vw::Vector2( llrModel[index].x(),
                                                                     llrModel[index].y() ) );
    xyzMatch[index] = vw::Vector3(0,0,0);
    float minDistance = std::numeric_limits<float>::max();
	vw::BBox2i region_to_search( DEM_pix - matchWindowHalfSize,
							 DEM_pix + matchWindowHalfSize + vw::Vector2i( 1, 1 ) );
	region_to_search.crop( DEM_bbox );
//	vw::Vector<int,2> start = DEM_pix - matchWindowHalfSize;
//    vw::Vector<int,2> limit = DEM_pix + matchWindowHalfSize + vw::Vector<int,2>( 1, 1 );
    
	// Get local copy from disk of the DEM that we are querying
	ImageView<float> DEM_window; // Assume this is small (less than 1024^2)
	if ( !IsMasked<typename ViewT::pixel_type>::value ) {
      DEM_window = crop(DEM.impl(), region_to_search );
	} else {
      DEM_window = crop(apply_mask(DEM.impl(), noDEMVal), region_to_search );
	}		

	// Ver2 of search in a neighborhood around x, y
	for( int k = 0; k < DEM_window.rows(); ++k ){
      for( int l = 0; l < DEM_window.cols(); ++l ){
        double interp = DEM_window( l,k );
        if( interp == noDEMVal ){ continue; }

        //compute the distance to the lidar point
        const vw::Vector2 lonlat = 
                DEMGeo.pixel_to_lonlat( Vector2(l,k) + region_to_search.min() );
        
        //transform into xyz coordinates of the foregound image.
        vw::Vector3 xyzDEM = DEMGeo.datum().geodetic_to_cartesian(  
            	                            vw::Vector3(lonlat.x(), lonlat.y(), interp) );
        xyzDEM = rotation*(xyzDEM-xyzMatchCenter) + xyzMatchCenter + translation;
        
        float distance = norm_2( xyzDEM - xyzModel[index] );

        //cout<<"distance="<<distance<<endl;

        if( distance < minDistance ){
          minDistance = distance;
          xyzMatch[index] = xyzDEM;
        }
      }
	}
    if( xyzMatch[index] == vw::Vector3(0,0,0) ){
      mask[index] = false;
    }
  }
  return mask;
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
                float matchError);



//used in DEM to DEM alignment
//old style DEM alignment
template <class ViewT>
void ICP_DEM_2_DEM(       std::vector<vw::Vector3>       featureArray, 
                    const vw::ImageViewBase<ViewT>&      backDEM,  
                    const vw::cartography::GeoReference& backDEMGeo, 
                    const vw::cartography::GeoReference& foreDEMGeo, 
		    const int matchingMode,
                    const CoregistrationParams&          settings,
                          vw::Vector3&                   translation, 
                          vw::Matrix<float, 3, 3>&       rotation, 
                          vw::Vector3&                   featureCenter, 
                          float&                         matchError ){

  std::vector<vw::Vector3> translationArray;
  std::vector<vw::Matrix<float, 3,3> > rotationArray;
  std::vector<vw::Vector3> referenceArray( featureArray.size() );
  vw::Vector3 referenceCenter;
  vw::Vector3 origFeatureCenter;
  vw::Vector3 zeroVector;

  std::vector<vw::Vector3> origFeatureArray( featureArray.size() );   
  for (unsigned int i = 0; i < featureArray.size(); i++){
    origFeatureArray[i] = featureArray[i];
  }
  origFeatureCenter = find_centroid( origFeatureArray );

  int numIter = 0;
  matchError = 100.0; 
   
  while( (numIter    < settings.maxNumIter) &&
         (matchError > settings.minConvThresh) ){
    std::cout << " iteration=" << numIter << std::endl;

    std::cout << "feature matching ..." << std::endl;
    std::vector<bool> mask = FindMatches( featureArray, backDEM, backDEMGeo, foreDEMGeo, 
                 referenceArray, settings.matchWindowHalfSize, settings.noDataVal);

    unsigned int num_good_matches = 0;
    for( std::vector<bool>::const_iterator i = mask.begin(); i != mask.end(); ++i ){
      if( *i ){ ++num_good_matches; }
    }

    std::vector<vw::Vector3> temp_feature;
    std::vector<vw::Vector3> temp_reference;
    if( num_good_matches == mask.size() ){
      temp_feature = featureArray;
      temp_reference = referenceArray;
    } else {
      temp_feature.reserve( num_good_matches );
      temp_reference.reserve( num_good_matches );
      for( unsigned int i=0; i < mask.size(); ++i ){
        if( mask[i] ){
          temp_feature.push_back( featureArray[i] );
          temp_reference.push_back( referenceArray[i] );
        }
      }
    }
    
    std::cout << "computing the matching error ..." << std::endl; 
    matchError = ComputeMatchingError3D( temp_feature, temp_reference );
    std::cout << "match error=" << matchError << std::endl;

    if (matchError >= std::numeric_limits<float>::max()){
      std::cout<< "Features cannot be matched. Exiting ICP_DEM_2_DEM..."<<std::endl;
      return;
    }

    featureCenter = find_centroid( temp_feature );
    referenceCenter = find_centroid( temp_reference );
    /*
    if (matchingMode >= 2){//ICP_TRANSLATION
      std::cout << "computing DEM translation ..." << std::endl;
      translation = ComputeDEMTranslation( featureArray, referenceArray );   
    } 
    
    std::cout << "translation = " << translation << std::endl;
    */
    
    if (matchingMode == 3){//ICP_ROTATION
      std::cout << "computing DEM rotation ..." << std::endl;
      //std::cout << "referenceCenter=" << referenceCenter << ", featureCenter=" << featureCenter << std::endl;
      rotation = ComputeDEMRotation( temp_feature, temp_reference, featureCenter, referenceCenter );
      
      //compute the rotation matrix with the origin of the referenceArray
      std::cout << "rotation matrix = " << std::endl;
      PrintMatrix( rotation );
    }
    
    if (matchingMode >= 2){//ICP_TRANSLATION
      std::cout << "computing DEM translation ..." << std::endl;
      //std::cout << "rotation matrix = " << std::endl;
      //PrintMatrix( rotation );
      //translation = ComputeDEMTranslation( temp_feature, temp_reference, rotation, featureCenter);  
      translation = ComputeDEMTranslation( temp_feature, temp_reference );  
      std::cout << "translation = " << translation << std::endl; 
    } 

    //apply the computed rotation and translation to the featureArray  
    //correction to transform features to ignore the center of rotation
    //TransformFeatures( featureArray, referenceCenter/*featureCenter*/, translation, rotation ); 

    
    //this is really bad
    Vector3 zeroVector;
    zeroVector(0)=0;
    zeroVector(1)=0;
    zeroVector(2)=0;
    //TransformFeatures( featureArray, zeroVector, translation, rotation ); 
    TransformFeatures( featureArray, featureCenter, translation, rotation ); 
    
    //save current rotation and translation
    translationArray.push_back( translation );
    rotationArray.push_back(    rotation    );
    
    ++numIter;
  }


  //compute the final rotation and translation - START
  rotation = rotationArray[0];
  translation = translationArray[0];
  for ( unsigned int i = 1; i < rotationArray.size(); ++i){
    rotation = rotationArray[i]*rotation;
    translation =  rotationArray[i]*translation + translationArray[i];
  }  
  //compute the final rotation and translation - END

  //for debug - START
  //cout<<"Alignment using the center of rotation"<<endl;
  //ICP_Report(origFeatureArray, featureArray, referenceArray, backDEMGeo, foreDEMGeo, featureCenter/*zeroVector*/, translation, rotation, matchError);
  //for debug - END
}


/*
//used in DEM to DEM alignment
//new style DEM alignment
template <class ViewT>
void ICP_DEM_2_DEM(       std::vector<vw::Vector3>       featureArray, 
                    const vw::ImageViewBase<ViewT>&      backDEM,  
                    const vw::cartography::GeoReference& backDEMGeo, 
                    const vw::cartography::GeoReference& foreDEMGeo, 
		    const int matchingMode,
                    const CoregistrationParams&          settings,
                          vw::Vector3&                   translation, 
                          vw::Matrix<float, 3, 3>&       rotation, 
                          vw::Vector3&                   featureCenter, 
                          float&                         matchError ){

  std::vector<vw::Vector3> translationArray;
  std::vector<vw::Matrix<float, 3,3> > rotationArray;
  std::vector<vw::Vector3> referenceArray( featureArray.size() );
  vw::Vector3 referenceCenter;
  vw::Vector3 origFeatureCenter;

  std::vector<vw::Vector3> origFeatureArray( featureArray.size() );   
  for (int i = 0; i < featureArray.size(); i++){
    origFeatureArray[i] = featureArray[i];
  }
  origFeatureCenter = find_centroid( origFeatureArray );

  int numIter = 0;
  matchError = 100.0; 
   
  while( (numIter    < settings.maxNumIter) &&
         (matchError > settings.minConvThresh) ){
    std::cout << " iteration=" << numIter << std::endl;

    std::cout << "feature matching ..." << std::endl;
    FindMatches( featureArray, backDEM, backDEMGeo, foreDEMGeo, 
                 referenceArray, settings.matchWindowHalfSize, settings.noDataVal);
    
    std::cout << "computing the matching error ..." << std::endl; 
    matchError = ComputeMatchingError3D( featureArray, referenceArray );
    std::cout << "match error=" << matchError << std::endl;

    if (matchError >= std::numeric_limits<float>::max()){
      std::cout<< "Features cannot be matched. Exiting ICP_DEM_2_DEM..."<<std::endl;
      return;
    }

    featureCenter = find_centroid( featureArray );
    referenceCenter = find_centroid(referenceArray);
    
    if (matchingMode >= 2){//ICP_TRANSLATION
      std::cout << "computing DEM translation ..." << std::endl;
      translation = ComputeDEMTranslation( featureArray, referenceArray );   
    } 
    std::cout << "translation = " << translation << std::endl;
    
    
    if (matchingMode == 3){//ICP_ROTATION
      std::cout << "computing DEM rotation ..." << std::endl;
      //std::cout << "referenceCenter=" << referenceCenter << ", featureCenter=" << featureCenter << std::endl;
      //rotation = ComputeDEMRotation( featureArray, referenceArray, featureCenter, referenceCenter );
      //compute the rotation matrix with the origin of the referenceArray
      rotation = ComputeDEMRotation( featureArray, referenceArray, translation);
    }
    std::cout << "rotation matrix = " << std::endl;
    PrintMatrix( rotation );

    Vector3 zeroVector;
    zeroVector(0)=0;
    zeroVector(1)=0;
    zeroVector(2)=0;

    //apply the computed rotation and translation to the featureArray  
    //correction to transform features to ignore the center of rotation
    TransformFeatures( featureArray, zeroVector, translation, rotation );     
        
    //save current rotation and translation
    translationArray.push_back( translation );
    rotationArray.push_back(    rotation    );
    
    ++numIter;
  }


  //compute the final rotation and translation - START
  rotation = rotationArray[0];
  translation = translationArray[0];
  for ( unsigned int i = 1; i < rotationArray.size(); ++i){
     rotation = rotationArray[i]*rotation;
     translation =  rotationArray[i]*translation + translationArray[i];
  }
  
  featureCenter(0)=0.0;
  featureCenter(1)=0.0;
  featureCenter(2)=0.0;
  
  cout<<"Alignment without center of rotation"<<endl;
  //compute the final rotation and translation - END
  ICP_Report(origFeatureArray, featureArray, referenceArray, backDEMGeo, foreDEMGeo, featureCenter, translation, rotation, matchError);
  
}
*/

//used in DEM to Lidar alignment
template <class ViewT>
void ICP_LIDAR_2_DEM(       std::vector<vw::Vector3>&      xyzMatchArray,  
                      const vw::ImageViewBase<ViewT>&      DEM,
                      const vw::cartography::GeoReference& DEMGeo, 
                      const std::vector<vw::Vector3>&      xyzModelArray, 
                      const std::vector<vw::Vector3>&      llrModelArray,
                      const CoregistrationParams           settings,
                            vw::Vector3&                   translation, 
                            vw::Matrix<float,3,3>&         rotation, 
                      const vw::Vector3&                   xyzModelCenter,
                            valarray<float>&               xyzErrorArray ){
  std::vector<vw::Vector3> translationArray;
  std::vector<vw::Matrix<float,3,3> > rotationArray;
  int numIter = 0;
  float avgMatchError = std::numeric_limits<float>::max();
  vw::Vector3 xyzMatchCenter = xyzModelCenter;

  while( (numIter       < settings.maxNumIter) && 
         (avgMatchError > settings.minConvThresh) ){
    vw_out(vw::InfoMessage, "icp") << " -- Iteration " << numIter << " --" << endl;
    
    std::vector<bool> mask = FindMatchesFromDEM( xyzModelArray, llrModelArray, DEM, DEMGeo, 
                        xyzMatchArray, translation, rotation, xyzMatchCenter, 
                        settings.noDataVal, settings.matchWindowHalfSize );

    unsigned int num_good_matches = 0;
    for( std::vector<bool>::const_iterator i = mask.begin(); i != mask.end(); ++i ){
      if( *i ){ ++num_good_matches; }
    }

    std::vector<vw::Vector3> temp_xyzMatch;
    std::vector<vw::Vector3> temp_xyzModel;
    if( num_good_matches == mask.size() ){
      temp_xyzMatch = xyzMatchArray;
      temp_xyzModel = xyzModelArray;
    } else {
      temp_xyzMatch.reserve( num_good_matches );
      temp_xyzModel.reserve( num_good_matches );
      for( unsigned int i=0; i < mask.size(); ++i ){
        if( mask[i] ){
          temp_xyzMatch.push_back( xyzMatchArray[i] );
          temp_xyzModel.push_back( xyzModelArray[i] );
        }
      }
    }
    
    vw_out(vw::InfoMessage, "icp") << "computing the matching error ... ";
    //xyzErrorArray = ComputeMatchingError(xyzMatchArray, xyzModelArray);
    //xyzErrorArray = ComputeMatchingError(temp_xyzModel, temp_xyzMatch);
    std::valarray<float> temp_xyzError( temp_xyzModel.size() );
    ComputeMatchingError( temp_xyzError, temp_xyzModel, temp_xyzMatch );
    avgMatchError = temp_xyzError.sum()/temp_xyzError.size();
    vw_out(vw::InfoMessage, "icp") << avgMatchError << endl;
    
    // Repopulate the xyzErrorArray from temp_xyzError
    if( temp_xyzError.size() == xyzErrorArray.size() ){
      xyzErrorArray = temp_xyzError;
    } else {
      for( unsigned int i=0, t=0; i < xyzErrorArray.size(); ++i ){
        if( mask[i] ){
          xyzErrorArray[i] = temp_xyzError[t];
          ++t;
        } else { xyzErrorArray[i] = std::numeric_limits<float>::max(); }
      }
    }


    //determine the xyzMatchCentroid
    vw::Vector3 xyzMatchCenter = find_centroid( temp_xyzMatch );

    vw_out(vw::InfoMessage, "icp") << "computing DEM translation ... ";
    translation = ComputeDEMTranslation(temp_xyzMatch, temp_xyzModel);
    vw_out(vw::InfoMessage, "icp") << translation << endl;

    vw_out(vw::InfoMessage, "icp") << "computing DEM rotation ... " << endl;;
    rotation = ComputeDEMRotation( temp_xyzMatch, temp_xyzModel, xyzMatchCenter, xyzModelCenter );
    vw_out(vw::InfoMessage, "icp") << rotation << endl;
    
    //apply the computed rotation and translation to the featureArray  
    //TransformFeatures(xyzMatchArray, translation, rotation);
    
    translationArray.push_back( translation );
    rotationArray.push_back( rotation );
    
    rotation = rotationArray[0];
    translation = translationArray[0];
    for( unsigned int i = 1; i < rotationArray.size(); ++i ){
      rotation = rotation*rotationArray[i];
      translation += translationArray[i];
    }
    ++numIter;
  }
}

#endif
