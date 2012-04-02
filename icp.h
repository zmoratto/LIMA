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


Vector2 back_2_fore_lonlat(Vector2 back_lon_lat);
Vector2 fore_2_back_lonlat(Vector2 fore_lon_lat);

//computes the translation between the foreground and background pixels
Vector3
ComputeDEMTranslation(
	const vector<Vector3>& featureArray, 
	const vector<Vector3>& matchArray);

void PrintMatrix(Matrix<float, 3, 3> &A);

//compute the matching error vector and overall average error 
valarray<float> 
ComputeMatchingError(
		const vector<Vector3>& featureArray, 
		const vector<Vector3>& matchArray);
float
ComputeMatchingError3D(const vector<Vector3>& modelArray, 
                       const vector<Vector3>& matchArray);
Matrix<float, 3, 3>
ComputeDEMRotation(
	const vector<Vector3>& featureArray, 
	const vector<Vector3>& matchArray,
	const Vector3&         matchCenter
	);

Matrix<float, 3, 3> 
ComputeDEMRotation(const vector<Vector3>& modelArray, 
		   const vector<Vector3>& matchArray,
		   const Vector3& modelCenter,
                   const Vector3& matchCenter);

inline
Matrix<float, 3, 3>
ComputeDEMRotation(
	const vector<Vector3>& featureArray, 
	const vector<Vector3>& matchArray)
{
const Vector3 matchCenter = find_centroid( matchArray );
return ComputeDEMRotation( featureArray, matchArray, matchCenter);
}


//applies a 3D rotation and transform to a DEM
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 translation, Matrix<float,3,3> rotation);
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 featureCenter, Vector3 translation, Matrix<float,3,3> rotation);

void ICP(vector<Vector3> featureArray, vector<Vector3> modelArray, CoregistrationParams settings,
	    Vector3 &translation, Matrix<float, 3, 3> &rotation, vector<float> &errorArray);

//compute a set of features from the fore image.
template <class ViewT>
vector<Vector3> 
GetFeatures(ImageViewBase<ViewT> const& foreImg, GeoReference const &foreGeo,
            ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo, 
            Vector2 samplingStep, Vector2 delta_lonlat, float noDataVal)
{
 
  cout<<"Get Features..."<<endl;

  int verStep = samplingStep(0);
  int horStep = samplingStep(1);
  vector<Vector3> featureArray;

  ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
  
  for (int j = 0; j < foreImg.impl().rows(); j=j+verStep){
     for (int i = 0; i < foreImg.impl().cols(); i=i+horStep){

       //determine the location of the corresponding background point for each valid fore point
       //if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != -noDataVal)  { 
       if ((isnan(foreImg.impl()(i,j))!=FP_NAN) && (foreImg.impl()(i,j)) != noDataVal)  {
	 Vector2 forePix(i,j);
	 Vector2 fore_lonlat = foreGeo.pixel_to_lonlat(forePix);
         fore_lonlat(0)=fore_lonlat(0)+delta_lonlat(0);
	 fore_lonlat(1)=fore_lonlat(1)+delta_lonlat(1);

         Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));
         Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat3);
      
         featureArray.push_back(fore_xyz);
       }
     }
  }

  cout<<"numFeatures="<<featureArray.size()<<endl;
  return featureArray;

};

//find the closest background points to the fore points.
template <class ViewT>
void
FindMatches(vector<Vector3> featureArray, ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo, 
            GeoReference const &foreGeo, vector<Vector3>& referenceArray, Vector2 matchWindowHalfSize, float noValData)
{
 
  cout<<" ***FIND MATCHES..."<<endl;

 int matchWindowHalfWidth = matchWindowHalfSize(0);
 int matchWindowHalfHeight = matchWindowHalfSize(1);

 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation()); 
 for (int i = 0; i < featureArray.size(); i++){

     //revert fore from cartesian to spherical coordinates
     Vector3 fore_lonlat3 = foreGeo.datum().cartesian_to_geodetic(featureArray[i]);
     Vector2 fore_lonlat;
 
     fore_lonlat(0) = fore_lonlat3(0);
     fore_lonlat(1) = fore_lonlat3(1);

     //Vector2 back_lonlat = fore_2_back_lonlat(fore_lonlat);
     Vector2 back_lonlat = fore_lonlat;

     Vector2 backPix = backGeo.lonlat_to_pixel(back_lonlat);	 
     float x = backPix(0);
     float y = backPix(1);

     Vector3 back_xyz;//invalid point
     float minDistance = 10000000000.0;

     for (int k =  y-matchWindowHalfHeight; k < y + matchWindowHalfHeight+1; k++){
        for (int l = x - matchWindowHalfWidth; l < x + matchWindowHalfWidth+1; l++){
 
          Vector2 backPix(l,k);
	  Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
           
          //revert to fore lon lat system
          Vector2 t_back_lonlat2 = back_lonlat;

          Vector3 t_back_lonlat3;
          t_back_lonlat3(0) = t_back_lonlat2(0);
          t_back_lonlat3(1) = t_back_lonlat2(1); 
          t_back_lonlat3(2) = interpBackImg(l,k);
     

          if (t_back_lonlat3(2) != noValData){
	    //transform into xyz coordinates of the foregound image.
	    Vector3 t_back_xyz= foreGeo.datum().geodetic_to_cartesian(t_back_lonlat3);             
      
	    float distance1 = t_back_xyz(0) - featureArray[i](0);
	    float distance2 = t_back_xyz(1) - featureArray[i](1);
	    float distance3 = t_back_xyz(2) - featureArray[i](2);
	    float distance = sqrt(distance1*distance1 + distance2*distance2 + distance3*distance3);   
	    if (distance < minDistance){
	      minDistance = distance;
              back_xyz = t_back_xyz;
	    }
            
	  }//endif (t_back_lonlat3(2) != noValData)
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
void
FindMatchesFromDEM(const vector<Vector3>&	xyzModelArray, 
		   const vector<Vector3>&	llrModelArray,
		   const ImageViewBase<ViewT>&	DEM, 
		   const GeoReference&		DEMGeo, 
		   vector<Vector3>&		xyzMatchArray, 
		   const Vector3&		translation, 
		   const Matrix<float,3,3>&	rotation, 
                   const Vector3&               xyzMatchCenter,
		   const double&		noDEMVal, 
		   const Vector2&		matchWindowHalfSize)
{
 ImageViewRef<float> interpDEM;
 if( IsMasked<typename ViewT::pixel_type>::value == 0 ){
   interpDEM = pixel_cast<float>(interpolate(edge_extend(DEM.impl(),ConstantEdgeExtension()),
					     BilinearInterpolation()) );
   
   //cout << "NOT masked" <<endl;
  }
 else{
   interpDEM = pixel_cast<float>(interpolate(edge_extend(apply_mask(DEM.impl()),ConstantEdgeExtension()),
					     BilinearInterpolation()) );
   //cout << "MASKED" <<endl;
 }
  
 //Vector3 xyzMatchCenter = find_centroid(xyzMatchArray);
 cout<<"XYZ_MATCH_CENTER="<<xyzMatchCenter<<endl;
 cout<<DEMGeo<<endl;  

 for (unsigned int index = 0; index < xyzModelArray.size(); index++){
   
   //const Vector3 lidar_lonlatrad = DEMGeo.datum().cartesian_to_geodetic(xyzModelArray[index]);
   //const Vector3 lidar_lonlatrad = xyz_to_lon_lat_radius(xyzModelArray[index] );
   //cout<<"lidar_alt="<<lidar_lonlatrad(2)<<endl;      
   //cout<<"llrModelArray="<<llrModelArray[index]<<endl;
   const Vector2 lidar_lonlat(llrModelArray[index](0), llrModelArray[index](1));
   const Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(lidar_lonlat);
   
   //const Vector2 DEM_pix = DEMGeo.lonlat_to_pixel( subvector(llrModelArray[index], 0, 2) );
   
   xyzMatchArray[index] = Vector3(0,0,0);
   float minDistance = 1000000.0;
   
  
   //search in a neigborhood around x,y and determine the best match to LOLA
   for(int k = DEM_pix.y() - matchWindowHalfSize.y(); k < DEM_pix.y() + matchWindowHalfSize.y() +1; k++){
     for(int l = DEM_pix.x() - matchWindowHalfSize.x(); l < DEM_pix.x() + matchWindowHalfSize.x() +1; l++){
       if ((l>=0) && (k>=0) && (l<DEM.impl().cols()) && (k<DEM.impl().rows())){
	 if (interpDEM(l,k)!=noDEMVal){
	   //compute the distance to the lidar point

	   const Vector2 pix(l,k);
  
	   const Vector2 lonlat = DEMGeo.pixel_to_lonlat( /*Vector2(l,k)*/pix );
	   
	   //revert to lon lat rad system
	   Vector3 lonlatrad (lonlat.x(),lonlat.y(),interpDEM(l,k)); // This z value is in meters, since DEMGeo.datum() is.
	   
	   //transform into xyz coordinates of the foregound image.
	   Vector3 xyzDEM = DEMGeo.datum().geodetic_to_cartesian(lonlatrad);
	   xyzDEM = rotation*(xyzDEM-xyzMatchCenter) + xyzMatchCenter + translation;
	   
	   Vector3 distance_vector = xyzDEM - xyzModelArray[index];
	   
	   float distance = norm_2( distance_vector );

	   //cout<<"distance="<<distance<<endl;

	   if (distance < minDistance){
	     minDistance = distance;
	     xyzMatchArray[index] = xyzDEM;
	   }
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
  vector<Vector3> xyz3DErrorArray;
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
    xyz3DErrorArray = ComputeMatchingError3D(xyzModelArray, xyzMatchArray);
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

 /*
//used in DEM to Lidar alignment
template <class ViewT>
void 
ICP_LIDAR_2_DEM(vector<Vector3>& 	    featureArray,  
		const ImageViewBase<ViewT>& DEM,
		const GeoReference& 	    DEMGeo, 
		const vector<Vector3>& 	    modelArray, 
		const vector<Vector3>& 	    modelArrayLatLon,
		const CoregistrationParams  settings,
		Vector3&		    translation, 
		Matrix<float, 3, 3>&        rotation, 
		const Vector3& 		    modelCenter,
		valarray<float>&	    errorArray) // probably don't need this anymore
{
     
  float matchError = std::numeric_limits<float>::max();
  cout<<"icp_lidar_2_dem:find matches"<<endl;     
  
  FindMatchesFromDEM(modelArray, modelArrayLatLon, DEM, DEMGeo, featureArray, 
		     translation, rotation, settings.noDataVal, settings.matchWindowHalfSize);
  
  vw_out(vw::InfoMessage, "icp") << "computing the matching error ... ";
  errorArray = ComputeMatchingError(featureArray, modelArray);
  
  
  cout<<"NUM_ERRORS="<<errorArray.size()<<endl;
  if (errorArray.size()!=0){
    matchError = errorArray.sum()/errorArray.size();
  }
  cout<<"matchError="<<matchError<<endl;
  vw_out(vw::InfoMessage, "icp") << matchError << endl;

  vector<Vector3> translationArray;
  vector<Matrix<float, 3,3> > rotationArray;
  //rotationArray.clear();
  //translationArray.clear();
  int numIter = 0;

  while((numIter < settings.maxNumIter) && (matchError > settings.minConvThresh)){
    vw_out(vw::InfoMessage, "icp") << " -- Iteration " << numIter << " --" << endl;
  
    vw_out(vw::InfoMessage, "icp") << "computing DEM translation ... ";
    translation = ComputeDEMTranslation(featureArray, modelArray);
    cout<<"translation="<<translation<<endl;
    vw_out(vw::InfoMessage, "icp") << translation << endl;
    
    vw_out(vw::InfoMessage, "icp") << "computing DEM rotation ... " << endl;;
    rotation = ComputeDEMRotation(featureArray, modelArray, modelCenter);
    vw_out(vw::InfoMessage, "icp") << rotation << endl;
    
    //apply the computed rotation and translation to the featureArray  
    TransformFeatures(featureArray, translation, rotation);
    
    translationArray.push_back(translation);
    rotationArray.push_back(rotation);
    
    rotation = rotationArray[0];
    translation = translationArray[0];
    for (unsigned int i = 1; i < rotationArray.size(); i++){
      rotation = rotation*rotationArray[i];
      translation += translationArray[i];
    }

    FindMatchesFromDEM(modelArray, modelArrayLatLon, DEM, DEMGeo, featureArray, 
		       translation, rotation, settings.noDataVal, 
		       settings.matchWindowHalfSize);
    
    vw_out(vw::InfoMessage, "icp") << "computing the matching error ... ";
    errorArray = ComputeMatchingError(featureArray, modelArray);
    matchError = errorArray.sum()/errorArray.size();
    cout<<"matchError="<<matchError<<endl;
    vw_out(vw::InfoMessage, "icp") << matchError << endl;
    
    numIter++;
  }
 
}
*/
  /*
#if 0
//determines best DEM points that match to LOLA track points
//used in ICP and in error calculations.
//the matching vector is returned in cartesian coordinates.
//the match is done by minimizing the Euclidian distance in the 3D cartesian space.
//TODO: add hor and ver pixel offset together with altitude
template <class ViewT>
void
FindMatchesFromDEM(const vector<Vector3>&	xyzModelArray, 
		   const vector<Vector3>&	llrModelArray,
		   const ImageViewBase<ViewT>&	DEM, 
		   const GeoReference&		DEMGeo, 
		   vector<Vector3>&		xyzMatchArray, 
                   const Vector3&               xyzMatchCenter,
		   const double&		noDEMVal, 
		   const Vector2&		matchWindowHalfSize)
{
 ImageViewRef<float> interpDEM;
 if( IsMasked<typename ViewT::pixel_type>::value == 0 ){
   interpDEM = pixel_cast<float>(interpolate(edge_extend(DEM.impl(),ConstantEdgeExtension()),
					     BilinearInterpolation()) );
   
   //cout << "NOT masked" <<endl;
  }
 else{
   interpDEM = pixel_cast<float>(interpolate(edge_extend(apply_mask(DEM.impl()),ConstantEdgeExtension()),
					     BilinearInterpolation()) );
   //cout << "MASKED" <<endl;
 }
  
 //Vector3 xyzMatchCenter = find_centroid(xyzMatchArray);
 cout<<"XYZ_MATCH_CENTER="<<xyzMatchCenter<<endl;
 cout<<DEMGeo<<endl;  

 for (unsigned int index = 0; index < xyzModelArray.size(); index++){
   
   //const Vector3 lidar_lonlatrad = DEMGeo.datum().cartesian_to_geodetic(xyzModelArray[index]);
   //const Vector3 lidar_lonlatrad = xyz_to_lon_lat_radius(xyzModelArray[index] );
   //cout<<"lidar_alt="<<lidar_lonlatrad(2)<<endl;      
   //cout<<"llrModelArray="<<llrModelArray[index]<<endl;
   const Vector2 lidar_lonlat(llrModelArray[index](0), llrModelArray[index](1));
   const Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(lidar_lonlat);
   
   //const Vector2 DEM_pix = DEMGeo.lonlat_to_pixel( subvector(llrModelArray[index], 0, 2) );
   
   xyzMatchArray[index] = Vector3(0,0,0);
   float minDistance = 1000000.0;
   
  
   //search in a neigborhood around x,y and determine the best match to LOLA
   for(int k = DEM_pix.y() - matchWindowHalfSize.y(); k < DEM_pix.y() + matchWindowHalfSize.y() +1; k++){
     for(int l = DEM_pix.x() - matchWindowHalfSize.x(); l < DEM_pix.x() + matchWindowHalfSize.x() +1; l++){
       if ((l>=0) && (k>=0) && (l<DEM.impl().cols()) && (k<DEM.impl().rows())){
	 if (interpDEM(l,k)!=noDEMVal){
	   //compute the distance to the lidar point

	   const Vector2 pix(l,k);
  
	   const Vector2 lonlat = DEMGeo.pixel_to_lonlat( pix );
	   
	   //revert to lon lat rad system
	   Vector3 lonlatrad (lonlat.x(),lonlat.y(),interpDEM(l,k)); // This z value is in meters, since DEMGeo.datum() is.
	   
	   //transform into xyz coordinates of the foregound image.
	   Vector3 xyzDEM = DEMGeo.datum().geodetic_to_cartesian(lonlatrad);
	   //xyzDEM = rotation*(xyzDEM-xyzMatchCenter) + xyzMatchCenter + translation;
	   
	   Vector3 distance_vector = xyzDEM - xyzModelArray[index];
	   
	   float distance = norm_2( distance_vector );

	   //cout<<"distance="<<distance<<endl;

	   if (distance < minDistance){
	     minDistance = distance;
	     xyzMatchArray[index] = xyzDEM;
	   }
	 }
       }
     }
   }
 }
}
#endif
*/

#if 0
template <class ViewT>
void
GetMatchesFromDEM(const vector<Vector3>&	xyzModelArray,
		  const vector<Vector3>&	llrModelArray,
		  const ImageViewBase<ViewT>&	DEM, 
		  const GeoReference&		DEMGeo, 
		  const double&		        noDEMVal,
		  vector<Vector3>&		xyzMatchArray,
		  valarray<float>&	        xyzErrorArray)
{

  ImageViewRef<int16> interpDEM;

  if( IsMasked<typename ViewT::pixel_type>::value == 0 ){
    interpDEM = pixel_cast<int16>(interpolate(edge_extend(DEM.impl(),ConstantEdgeExtension()),
					      BilinearInterpolation()) );    
    cout << "NOT masked" <<endl;
  }
  else{
    interpDEM = pixel_cast<int16>(interpolate(edge_extend(apply_mask(DEM.impl()),ConstantEdgeExtension()),
					      BilinearInterpolation()) );
    cout << "MASKED" <<endl;
  }
  
  cout<<"GetMatchesFromDEM: numLOLAPoints="<<llrModelArray.size()<<endl;
  
  for (unsigned int index = 0; index < llrModelArray.size(); index++){
    
   //const Vector3 lidar_lonlatrad = DEMGeo.datum().cartesian_to_geodetic(lidar_xyz[index]);
   //const Vector3 lidar_lonlatrad = xyz_to_lon_lat_radius( lidar_xyz[index] );
   //cout<<"lidar_alt="<<lidar_lonlatrad(2)<<endl;      
   //const Vector2 lidar_lonlat(lidar_lonlatrad(0), lidar_lonlatrad(1));
   //cout<<"test1"<<endl;
    
    const Vector2 DEM_pix = DEMGeo.lonlat_to_pixel( subvector(llrModelArray[index], 0, 2) );
    //cout<<DEM_pix<<", width="<<interpDEM.cols()<<", height="<<interpDEM.rows()<<endl;

    //set the default values
    xyzMatchArray[index] = Vector3(0,0,0);
    xyzErrorArray[index] = -1;
    
    float l = DEM_pix.x(); 
    float k = DEM_pix.y(); 
    //cout<<"l="<<l<<", k="<<k<<endl;

    if ((l>=0) && (k>=0) && (l<interpDEM.cols()) && (k<interpDEM.rows())){
      //if (interpDEM(l,k) != noDEMVal){
      if ((interpDEM(l,k) < 4000) && (interpDEM(l,k) > -4000)){
	/*
	//compute the distance to the lidar point
	//const Vector2 pix(l,k);
	const Vector2 lonlat = DEMGeo.pixel_to_lonlat( Vector2(l,k) );

	//revert to lon lat rad system
        
	Vector3 lonlatrad (lonlat.x(),lonlat.y(),interpDEM(l,k)); // This z value is in meters, since DEMGeo.datum() is.
	*/
        Vector3 lonlatrad (llrModelArray[index](0), llrModelArray[index](1), interpDEM(l,k)); // This z value is in meters, since DEMGeo.datum() is.

	//transform into xyz coordinates of the foregound image.
	Vector3 dem_xyz = DEMGeo.datum().geodetic_to_cartesian(lonlatrad);
       
	Vector3 distance_vector = dem_xyz - xyzModelArray[index];
	
	float distance = norm_2( distance_vector );
	
	xyzMatchArray[index] = dem_xyz;
	xyzErrorArray[index] = distance;
	
	//cout<<"GetMatchesFromDEM: lidar="<<llrModelArray[index][2]<<", dem="<<1737400+interpDEM(l,k)<<", error="<<radErrorArray[index]<<endl;
	
      }
    }
  }
  
}
#endif

#endif
