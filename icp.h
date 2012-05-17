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

void PrintMatrix(Matrix<float, 3, 3> &A);

//compute the matching error vector and overall average error 
valarray<float> ComputeMatchingError( const std::vector<vw::Vector3>& model, 
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
template <class ViewT>
void FindMatches( const std::vector<vw::Vector3>&      features, 
                  const vw::ImageViewBase<ViewT>&      backImg, 
                  const vw::cartography::GeoReference& backGeo, 
                  const vw::cartography::GeoReference& foreGeo, 
                        std::vector<vw::Vector3>&      referenceArray, 
                  const vw::Vector2&                   matchWindowHalfSize, 
                  const float                          noValData) {
  //cout<<" ***FIND MATCHES..."<<endl;
  vw::ImageViewRef<typename ViewT::pixel_type> interpBackImg = 
        interpolate(edge_extend(backImg.impl(), vw::ConstantEdgeExtension()), 
                    vw::BilinearInterpolation()); 

  for( unsigned int i = 0; i < features.size(); ++i ){
    //revert fore from cartesian to spherical coordinates
    vw::Vector3 fore_lonlat3 = foreGeo.datum().cartesian_to_geodetic( features[i] );

    vw::Vector2 backPix = backGeo.lonlat_to_pixel( vw::Vector2( fore_lonlat3.x(), 
                                                                fore_lonlat3.y() ) );

    vw::Vector3 back_xyz; //invalid point
    float minDistance = std::numeric_limits<float>::max();
    vw::Vector<int,2> start = backPix - matchWindowHalfSize;
    vw::Vector<int,2> limit = backPix + matchWindowHalfSize + vw::Vector<int,2>( 1, 1 );

    for( int k = start.y(); k < limit.y(); ++k ){
      for( int l = start.x(); l < limit.x(); ++l ){
        float interpBack = interpBackImg(l,k);
        if( interpBack == noValData ){ continue; }

        vw::Vector2 bPixel = backGeo.pixel_to_lonlat( vw::Vector2(l,k) );
          
        //revert to fore lon lat system
        vw::Vector3 t_back_lonlat3( bPixel.x(), bPixel.y(), interpBack );

        //transform into xyz coordinates of the foregound image.
        vw::Vector3 t_back_xyz= foreGeo.datum().geodetic_to_cartesian(t_back_lonlat3);
        float distance = norm_2( t_back_xyz - features[i] );
     
        if( distance < minDistance ){
          minDistance = distance;
          back_xyz = t_back_xyz;
        }
      }
    }
    referenceArray[i]=back_xyz;
  }
};

//determines best DEM points that match to LOLA track points
//used in ICP and in error calculations.
//the matching vector is returned in cartesian coordinates.
//the match is done by minimizing the Euclidian distance in the 3D cartesian space.
//TODO: add hor and ver pixel offset together with altitude
template <class ViewT>
void FindMatchesFromDEM( const std::vector<vw::Vector3>&      xyzModel, 
                         const std::vector<vw::Vector3>&      llrModel,
                         const vw::ImageViewBase<ViewT>&      DEM, 
                         const vw::cartography::GeoReference& DEMGeo, 
                               std::vector<vw::Vector3>&      xyzMatch, 
                         const vw::Vector3&                   translation, 
                         const vw::Matrix<float,3,3>&         rotation, 
                         const vw::Vector3&                   xyzMatchCenter,
                         const double                         noDEMVal, 
                         const vw::Vector2&                   matchWindowHalfSize ) {
  vw::ImageViewRef<float> interpDEM;
  if( IsMasked<typename ViewT::pixel_type>::value == 0 ){
    interpDEM = pixel_cast<float>(interpolate( edge_extend(DEM.impl(),
                                                           vw::ConstantEdgeExtension()),
                                               vw::BilinearInterpolation()) );
    //cout << "NOT masked" <<endl;
  }
  else{
    interpDEM = pixel_cast<float>(interpolate( edge_extend( apply_mask(DEM.impl()),
                                                            vw::ConstantEdgeExtension()),
                                               vw::BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }
  
  // cout<<"XYZ_MATCH_CENTER="<<xyzMatchCenter<<endl;
  // cout<<DEMGeo<<endl;  

  vw::BBox2 DEM_bbox = bounding_box( DEM );
  for( unsigned int index = 0; index < xyzModel.size(); ++index ){
    const vw::Vector2 DEM_pix = DEMGeo.lonlat_to_pixel( vw::Vector2( llrModel[index].x(),
                                                                     llrModel[index].y() ) );
    xyzMatch[index] = vw::Vector3(0,0,0);
    float minDistance = std::numeric_limits<float>::max();
    vw::Vector<int,2> start = DEM_pix - matchWindowHalfSize;
    vw::Vector<int,2> limit = DEM_pix + matchWindowHalfSize + vw::Vector<int,2>( 1, 1 );
    
    //search in a neigborhood around x,y and determine the best match to LOLA
    for( int k = start.y(); k < limit.y(); ++k ){
      for( int l = start.x(); l < limit.x(); ++l){
        if( DEM_bbox.contains(vw::Vector<int,2>( l,k )) ){
          double interp = interpDEM(l,k);
          if( interp == noDEMVal ){ continue; }

          //compute the distance to the lidar point
          const vw::Vector2 lonlat = DEMGeo.pixel_to_lonlat( Vector2(l,k) );
        
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
    }
  }
}



//used in DEM to DEM alignment
//TODO: extract the xyz centroid position
template <class ViewT>
void 
ICP_DEM_2_DEM(vector<Vector3> featureArray, ImageViewBase<ViewT> const& backDEM,  
              GeoReference const& backDEMGeo, GeoReference const& foreDEMGeo, CoregistrationParams settings,
              Vector3 &translation, Matrix<float, 3, 3> &rotation, Vector3 &featureCenter, /*vector<float> &errorArray*/float &matchError )
{
   
    vector<Vector3> translationArray;
    vector<Matrix<float, 3,3> > rotationArray;
    vector<Vector3> referenceArray;
    referenceArray.resize(featureArray.size());
    
    int numIter = 0;
    matchError = 100.0; 
   
    while((numIter < settings.maxNumIter)&&(matchError > settings.minConvThresh)){
            
            cout<<"iteration="<<numIter<<endl;

	    printf("feature matching ...\n");
	    FindMatches(featureArray, backDEM, backDEMGeo, foreDEMGeo, referenceArray, 
                        settings.matchWindowHalfSize, settings.noDataVal);
      
	    cout<<"computing the matching error ..."<<endl; 
            matchError = ComputeMatchingError3D(featureArray, referenceArray);
            cout<<"match error="<<matchError<<endl;

	    featureCenter = find_centroid(featureArray);
	    Vector3 referenceCenter = find_centroid(referenceArray);

	    cout<<"computing DEM translation ..."<<endl;
	    translation = ComputeDEMTranslation(featureArray, referenceArray);
            cout<<"translation = "<<translation<<endl;
             
	    cout<<"computing DEM rotation ..."<<endl;
	    cout<<"referenceCenter="<<referenceCenter<<", featureCenter="<<featureCenter<<endl;
	    
	    rotation = ComputeDEMRotation(featureArray, referenceArray, featureCenter, referenceCenter);
	    cout<<"rotation matrix = "<<endl;
	    PrintMatrix(rotation);

	    //apply the computed rotation and translation to the featureArray  
	    TransformFeatures(featureArray, featureCenter, translation, rotation);
            
            //save current rotation and translation
	    translationArray.push_back(translation);
	    rotationArray.push_back(rotation);
	    
	    numIter++;
           
    }
    //compute the final rotation and translation
    rotation = rotationArray[0];
    translation = translationArray[0];
    for (int i = 1; i < rotationArray.size(); i++){
	 rotation = rotationArray[i]*rotation;
	 translation = translationArray[i] + rotationArray[i]*translation;
    }
 
}


//used in DEM to Lidar alignment
template <class ViewT>
void 
ICP_LIDAR_2_DEM(vector<Vector3>& 	    xyzMatchArray,  
		const ImageViewBase<ViewT>& DEM,
		const GeoReference& 	    DEMGeo, 
		const vector<Vector3>& 	    xyzModelArray, 
		const vector<Vector3>& 	    llrModelArray,
		const CoregistrationParams  settings,
		Vector3&		    translation, 
		Matrix<float, 3, 3>&        rotation, 
		const Vector3& 		    xyzModelCenter,
		valarray<float>&	    xyzErrorArray)
{
   
  vector<Vector3> translationArray;
  vector<Matrix<float, 3,3> > rotationArray;
  //vector<Vector3> xyz3DErrorArray;
  int numIter = 0;
  float avgMatchError = std::numeric_limits<float>::max();
  //rotationArray.clear();
  //translationArray.clear();
  Vector3 xyzMatchCenter = xyzModelCenter;

  while((numIter < settings.maxNumIter) && (avgMatchError > settings.minConvThresh)){
    vw_out(vw::InfoMessage, "icp") << " -- Iteration " << numIter << " --" << endl;
    
    FindMatchesFromDEM(xyzModelArray, llrModelArray, DEM, DEMGeo, xyzMatchArray, 
		       translation, rotation, xyzMatchCenter, settings.noDataVal, 
		       settings.matchWindowHalfSize);
    
    vw_out(vw::InfoMessage, "icp") << "computing the matching error ... ";
    //xyz3DErrorArray = ComputeMatchingError3D(xyzModelArray, xyzMatchArray);
    xyzErrorArray = ComputeMatchingError(xyzMatchArray, xyzModelArray);
    avgMatchError = xyzErrorArray.sum()/xyzErrorArray.size();
    cout<<"avgMatchError="<<avgMatchError<<endl;
    vw_out(vw::InfoMessage, "icp") << avgMatchError << endl;

    //determine the xyzMatchCentroid
    Vector3 xyzMatchCenter = find_centroid(xyzMatchArray);

    vw_out(vw::InfoMessage, "icp") << "computing DEM translation ... ";
    translation = ComputeDEMTranslation(xyzMatchArray, xyzModelArray);
    cout<<"translation="<<translation<<endl;
    vw_out(vw::InfoMessage, "icp") << translation << endl;

    vw_out(vw::InfoMessage, "icp") << "computing DEM rotation ... " << endl;;
    rotation = ComputeDEMRotation(xyzMatchArray, xyzModelArray, /*xyzModelCenter*/xyzMatchCenter, xyzModelCenter);
    vw_out(vw::InfoMessage, "icp") << rotation << endl;
    
    //apply the computed rotation and translation to the featureArray  
    //TransformFeatures(xyzMatchArray, translation, rotation);
    
    translationArray.push_back(translation);
    rotationArray.push_back(rotation);
    
    rotation = rotationArray[0];
    translation = translationArray[0];
    for (unsigned int i = 1; i < rotationArray.size(); i++){
      rotation = rotation*rotationArray[i];
      translation += translationArray[i];
    }
    
    numIter++;
  }
 
}

#endif
