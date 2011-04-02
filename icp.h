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


//computes the translation between the foreground and background pixels
void ComputeDEMTranslation(vector<Vector3> featureArray, vector<Vector3> matchArray, Vector3 &translation );

void PrintMatrix(Matrix<float, 3, 3> &A);

//compute the matching error vector and overall average error 
float ComputeMatchingError(vector<Vector3> featureArray, vector<Vector3> matchArray, vector<float> &errorArray);

void ComputeDEMRotation(vector<Vector3> featureArray, vector<Vector3> matchArray, 
                        Vector3 translation, Matrix<float, 3, 3> &rotation);

//applies a 3D rotation and transform to a DEM
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 translation, Matrix<float,3,3> rotation);




//compute a set of features from the fore image.
template <class ViewT>
vector<Vector3> 
GetFeatures(ImageViewBase<ViewT> const& foreImg, GeoReference const &foreGeo,
            ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo)
{
 
 int verStep = 8;
 int horStep = 8;
 vector<Vector3> featureArray;
 //determine if this is ppd or mpp.
 //determine the modata value as well
 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
  
 for (int j = 0; j < foreImg.impl().rows(); j=j+verStep){
     for (int i = 0; i < foreImg.impl().cols(); i=i+horStep){
       //determine the location of the corresponding background point

       if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != 0)  { 

	 Vector2 forePix(i,j);
	 Vector2 fore_lonlat = foreGeo.pixel_to_lonlat(forePix);
	 //fore_lonlat(0) = fore_lonlat(0)-180;
       
         Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));
         Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat3);
      
         featureArray.push_back(fore_xyz);
       }
     }
  }

 return featureArray;

};

//find the closest background points to the fore points.
template <class ViewT>
void
FindMatches(vector<Vector3> featureArray, ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo, 
            GeoReference const &foreGeo, vector<Vector3>& matchArray )
{
 
 int matchWindowHalfWidth = 5;
 int matchWindowHalfHeight = 5;
 float usgs_2_lonlat = 180/(3.14159265*3396190);

 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 
 for (int i = 0; i < featureArray.size(); i++){

     //TO DO:revert fore from cartesian to spherical coordinates
     Vector3 fore_lonlat3 = foreGeo.datum().cartesian_to_geodetic(featureArray[i]);
     Vector2 fore_lonlat;
     fore_lonlat(0) = fore_lonlat3(0)-180; 
     fore_lonlat(1) = fore_lonlat3(1);
     
     Vector2 back_lonlat = fore_lonlat/usgs_2_lonlat;
     Vector2 backPix = backGeo.lonlat_to_pixel(back_lonlat);	 
     float x = backPix(0);
     float y = backPix(1);
     
     //initialize the search with the closest position on the background DEM
     //Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
     Vector3 back_lonlat3;
     back_lonlat3(0) = back_lonlat(0)*usgs_2_lonlat+180;
     back_lonlat3(1) = back_lonlat(1)*usgs_2_lonlat; 
     back_lonlat3(2) = interpBackImg(x,y);
     
     //cout<<"fore: "<<fore_lonlat3<<endl;
     //cout<<"back_direct: "<<back_lonlat3<<endl;
     

     Vector3 back_xyz;
     //this can be commented out - START
     float minDistance = 10000000000.0;
     for (int k =  y-matchWindowHalfHeight; k < y + matchWindowHalfHeight; k++){
        for (int l = x - matchWindowHalfWidth; l < x + matchWindowHalfWidth; l++){
 
          Vector2 backPix(l,k);
	  Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
          Vector3 t_back_lonlat3;
          
          //revert to fore lon lat system
          //t_back_lonlat3(0) = back_lonlat(0)*usgs_2_lonlat; 
          t_back_lonlat3(0) = back_lonlat(0)*usgs_2_lonlat+180;
          t_back_lonlat3(1) = back_lonlat(1)*usgs_2_lonlat; 
          t_back_lonlat3(2) = interpBackImg(l,k);

          //cout<<backPix<<endl;
          //cout<<"fore: "<<fore_lonlat3<<endl;          
          //cout<<"back: "<<t_back_lonlat3<<endl;

          //transform into xyz coordinates of the foregound image.
          Vector3 t_back_xyz= foreGeo.datum().geodetic_to_cartesian(t_back_lonlat3);             
          //cout<<"fore"<<featureArray[i](0)<<" "<<featureArray[i](1)<<" "<<featureArray[i](2)<<endl;          
          //cout<<"back"<<t_back_lonlat3(0)<<" "<<t_back_lonlat3(1)<<" "<<t_back_lonlat3(2)<<endl;

          float distance1 = t_back_xyz(0) - featureArray[i](0);
          float distance2 = t_back_xyz(1) - featureArray[i](1);
          float distance3 = t_back_xyz(2) - featureArray[i](2);
          float distance = sqrt(distance1*distance1 + distance2*distance2 + distance3*distance3);   
          if (distance < minDistance){
	    minDistance = distance;
            back_xyz = t_back_xyz;
          }
  
	}
     }
     //this can be commented out - END
     //cout<<"back_window: "<<back_lonlat3<<endl;
    
     matchArray[i]=back_xyz;
 }
 
};
