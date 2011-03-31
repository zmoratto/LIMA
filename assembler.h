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

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
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




template <class ViewT1, class ViewT2 >
void
ComputeAssembledImage(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                      ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                      string assembledImgFilename, int mode, Vector3 translation, Matrix<float,3,3> rotation)
{
 
  float usgs_2_lonlat = 180/(3.14159265*3396190);
  float maxUpsampleRatioBackImg = 4.0;

  float numPixPerDegreeBackImg = ComputePixPerDegree(backGeo, orig_backImg.impl().cols(), orig_backImg.impl().rows(), 1);
  printf("numPixPerDegreeBackmg = %f\n", numPixPerDegreeBackImg);
  float numPixPerDegreeForeImg = ComputePixPerDegree(foreGeo, orig_foreImg.impl().cols(), orig_foreImg.impl().rows(), 0);
  printf("numPixPerDegreeForeImg = %f\n", numPixPerDegreeForeImg);
  
  float upsampleRatioBackImg = numPixPerDegreeForeImg/numPixPerDegreeBackImg;
  float upsampleRatioForeImg = 1;
  if (upsampleRatioBackImg > maxUpsampleRatioBackImg){
      upsampleRatioForeImg = maxUpsampleRatioBackImg/upsampleRatioBackImg;
      upsampleRatioBackImg = maxUpsampleRatioBackImg;  
  }
  
  printf("upsampleRatioBackDEM = %f\n", upsampleRatioBackImg);
  printf("upsampleRatioForeDEM = %f\n", upsampleRatioForeImg);

  int assembledWidth = (int)floor(upsampleRatioBackImg*orig_backImg.impl().cols());
  int assembledHeight = (int)floor(upsampleRatioBackImg*orig_backImg.impl().rows());
  printf("W = %d, H = %d, aW = %d, aH = %d\n", orig_foreImg.impl().cols(), orig_foreImg.impl().rows(), assembledWidth, assembledHeight);
  

  Matrix<double> H;
  ImageView<typename ViewT1::pixel_type> assembledImg(assembledWidth, assembledHeight);

  ImageViewRef<typename ViewT1::pixel_type>  foreImg = interpolate(edge_extend(orig_foreImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
   
  ImageViewRef<typename ViewT2::pixel_type>  backImg = interpolate(edge_extend(orig_backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 
  H = backGeo.transform();
  H(0,0) /=upsampleRatioBackImg;
  H(1,1) /=upsampleRatioBackImg;
  GeoReference assembledGeo; 
  assembledGeo.set_transform(H);

  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
  cout << "num_channels = "<< num_channels << "\n";

  for (int j = 0; j < assembledHeight-1; j++){
    for (int i = 0; i < assembledWidth-1; i++){
      
       float y = j/upsampleRatioBackImg;
       float x = i/upsampleRatioBackImg;     

       //copy the background image
       if (mode == 0){ //DEM
          assembledImg.impl()(i,j)[0] = backImg.impl()(x,y)[0];
       }
       else{ //DRG
	 assembledImg.impl()(i,j)[0] = backImg.impl()(x,y)[0];
	 assembledImg.impl()(i,j)[1] = backImg.impl()(x,y)[0];
	 assembledImg.impl()(i,j)[2] = backImg.impl()(x,y)[0];
       }

       //get the pixel with these coordinates in the foreDRG 
       Vector2 pix(x, y);
       Vector2 backUSGS_lonlat = backGeo.pixel_to_lonlat(pix);
       Vector2 back_lonlat = backUSGS_lonlat*usgs_2_lonlat;
           
       if (back_lonlat(0)< 0) back_lonlat(0) = 180+back_lonlat(0);
      
       Vector2 fore_pix = foreGeo.lonlat_to_pixel(back_lonlat);
       x = fore_pix(0);
       y = fore_pix(1);
       
       //copy the foreground image
       if ((x>0) && (x<foreImg.impl().cols()-1) && (y>0) && (y<foreImg.impl().rows())){

         if (mode == 0){
	   if ((isnan(foreImg.impl()(x,y)[0])!=FP_NAN) && (foreImg.impl()(x,y)[0]) != 0)  { 
             assembledImg.impl()(i,j) = foreImg.impl()(x,y) - translation[2]; 
	   }
	 }
         else{
	   if ((foreImg.impl()(x,y)[0]!=255) && (foreImg.impl()(x,y)[1]!=255) && (foreImg.impl()(x,y)[2]!=255)){
             assembledImg.impl()(i,j)[0] = foreImg.impl()(x,y)[0];
             assembledImg.impl()(i,j)[1] = foreImg.impl()(x,y)[1];
             assembledImg.impl()(i,j)[2] = foreImg.impl()(x,y)[2];
	   }
	 }
       }
       
   }
  }
  
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}



template <class ViewT1, class ViewT2 >
void
ComputeAssembledImageNew(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                         ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                         string assembledImgFilename, int mode, Vector3 translation, Matrix<float,3,3> rotation)
{
 
  cout<<"NEW STYLE"<<endl;

  float usgs_2_lonlat = 180/(3.14159265*3396190);
  float maxUpsampleRatioBackImg = 4.0;

  float numPixPerDegreeBackImg = ComputePixPerDegree(backGeo, orig_backImg.impl().cols(), orig_backImg.impl().rows(), 1);
  //printf("numPixPerDegreeBackmg = %f\n", numPixPerDegreeBackImg);
  float numPixPerDegreeForeImg = ComputePixPerDegree(foreGeo, orig_foreImg.impl().cols(), orig_foreImg.impl().rows(), 0);
  //printf("numPixPerDegreeForeImg = %f\n", numPixPerDegreeForeImg);
  
  float upsampleRatioBackImg = numPixPerDegreeForeImg/numPixPerDegreeBackImg;
  float upsampleRatioForeImg = 1;
  if (upsampleRatioBackImg > maxUpsampleRatioBackImg){
      upsampleRatioForeImg = maxUpsampleRatioBackImg/upsampleRatioBackImg;
      upsampleRatioBackImg = maxUpsampleRatioBackImg;  
  }
  
  printf("upsampleRatioBackDEM = %f\n", upsampleRatioBackImg);
  printf("upsampleRatioForeDEM = %f\n", upsampleRatioForeImg);

  int assembledWidth = (int)floor(upsampleRatioBackImg*orig_backImg.impl().cols());
  int assembledHeight = (int)floor(upsampleRatioBackImg*orig_backImg.impl().rows());
  printf("W = %d, H = %d, aW = %d, aH = %d\n", orig_foreImg.impl().cols(), orig_foreImg.impl().rows(), assembledWidth, assembledHeight);
  

  Matrix<double> H;
  ImageView<typename ViewT1::pixel_type> assembledImg(assembledWidth, assembledHeight);

  ImageViewRef<typename ViewT1::pixel_type>  foreImg = interpolate(edge_extend(orig_foreImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
   
  ImageViewRef<typename ViewT2::pixel_type>  backImg = interpolate(edge_extend(orig_backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 
  H = backGeo.transform();
  H(0,0) /=upsampleRatioBackImg;
  H(1,1) /=upsampleRatioBackImg;
  GeoReference assembledGeo; 
  assembledGeo.set_transform(H);

  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
  cout << "num_channels = "<< num_channels << "\n";

  //copy the background image to the assembled image.

  for (int j = 0; j < assembledHeight-1; j++){
    for (int i = 0; i < assembledWidth-1; i++){
      
       float y = j/upsampleRatioBackImg;
       float x = i/upsampleRatioBackImg;     

       //copy the background image
       if (mode == 0){ //DEM
          assembledImg.impl()(i,j)[0] = backImg.impl()(x,y)[0];
       }
       else{ //DRG
	 assembledImg.impl()(i,j)[0] = backImg.impl()(x,y)[0];
	 assembledImg.impl()(i,j)[1] = backImg.impl()(x,y)[0];
	 assembledImg.impl()(i,j)[2] = backImg.impl()(x,y)[0];
       }
    }
   }

  //transform the foreground image and copy it over the background image
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){
       
        if (mode == 0){
	  if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != 0)  {
             
             //get the 3D coordinates
	     Vector2 forePix(i,j);
             Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
             Vector3 fore_lon_lat_rad;
             fore_lon_lat_rad(0) = fore_lon_lat(0);
	     fore_lon_lat_rad(1) = fore_lon_lat(1);
             fore_lon_lat_rad(2) = foreImg.impl()(i,j)[0];

             Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lon_lat_rad);

             //transform the coordinates using the rotation/translation params
             Vector3 transf_fore_xyz = rotation*fore_xyz + translation;  
	      
             //transform in spherical coords
             Vector3 transf_fore_lon_lat_rad = foreGeo.datum().cartesian_to_geodetic(transf_fore_xyz);

             //TO DO:determine their location on the assembled image
             //TO DO:overwrite the assembled image 
	     //assembledImg.impl()(i,j) = foreImg.impl()(x,y) - translation[2];
	   }
        }
        else{
	  if ((foreImg.impl()(i,j)[0]!=255) && (foreImg.impl()(i,j)[1]!=255) && (foreImg.impl()(i,j)[2]!=255)){
              //get the coordinates
              Vector2 forePix(i,j);
              Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);

              //change into USGS coords
              fore_lon_lat(0) = fore_lon_lat(0)-180;
              //cout<<"fore init:lon="<<fore_lon_lat(0)<<",lat="<<fore_lon_lat(1)<<endl;
              fore_lon_lat = fore_lon_lat/usgs_2_lonlat;  
	      //cout<<"fore:lon="<<fore_lon_lat(0)<<",lat="<<fore_lon_lat(1)<<endl;
              
              //determine their location on the assembled image
              Vector2 backPix= assembledGeo.lonlat_to_pixel(fore_lon_lat);
              int x = backPix(0);
              int y = backPix(1);
              //cout<<"x="<<x<<", "<<"y="<<y<<",width="<<assembledImg.impl().cols()<<",height="<<assembledImg.impl().rows()<<endl;
              if ((x>0) && (y>0) && (x<assembledImg.impl().cols()) && (y<assembledImg.impl().rows())){ 
		//overwrite the assembled image 
		assembledImg.impl()(x,y)[0] = foreImg.impl()(i,j)[0];
		assembledImg.impl()(x,y)[1] = foreImg.impl()(i,j)[1];
		assembledImg.impl()(x,y)[2] = foreImg.impl()(i,j)[2];
	      }
	   }
	}
    }
  }

  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}


