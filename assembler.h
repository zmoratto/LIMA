// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include "icp.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;


//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDEM(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    string assembledImgFilename, Vector3 translation, Matrix<float,3,3>rotation, 
                    Vector3 center, Vector2 bestDeltaLonLat)
{
 
  float maxUpsampleRatioBackImg = 40.0;
  int maxAssembledImgWidth = 4096;
  int maxAssembledImgHeight = 4096;

  cout<<"Background"<<endl;
  float numPixPerDegreeBackImg = ComputePixPerDegree(backGeo, orig_backImg.impl().cols(), orig_backImg.impl().rows(), 1);
  //printf("numPixPerDegreeBackmg = %f\n", numPixPerDegreeBackImg);
   cout<<"Foreground"<<endl;
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

  if (assembledWidth > maxAssembledImgWidth){
      assembledWidth = maxAssembledImgWidth;
  }
  if (assembledHeight > maxAssembledImgHeight){
      assembledHeight = maxAssembledImgHeight;
  }
  printf("W = %d, H = %d, aW = %d, aH = %d\n", 
         orig_foreImg.impl().cols(), orig_foreImg.impl().rows(), assembledWidth, assembledHeight);
  
  Matrix<double> H;
  ImageView<typename ViewT2::pixel_type> assembledImg(assembledWidth, assembledHeight);
 
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
    = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
 
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
    = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  //compute the foreground centroid
  Vector3 foreCenter_xyz;
  Vector3 foreCenter_lon_lat_rad;
  int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){ 
        if ((isnan(foreImg.impl()(i,j))!=FP_NAN) && (foreImg.impl()(i,j)) != 0)  {       
            
	  //case of cartesian coordinates
	  //get the coordinates
	  Vector2 forePix(i,j);
	  Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	  fore_lon_lat(0)=fore_lon_lat(0)+bestDeltaLonLat(0);
	  fore_lon_lat(1)=fore_lon_lat(1)+bestDeltaLonLat(1);
	  
	  //spherical to cartesian
	  Vector3 fore_lon_lat_rad;
	  fore_lon_lat_rad(0) = fore_lon_lat(0);
	  fore_lon_lat_rad(1) = fore_lon_lat(1);
	  fore_lon_lat_rad(2) = foreImg.impl()(i,j);
	  Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lon_lat_rad); 
	  foreCenter_xyz = foreCenter_xyz + fore_xyz;
	  foreCenter_lon_lat_rad = foreCenter_lon_lat_rad + fore_lon_lat_rad;
	  count++;
	}
    }
  }

  foreCenter_xyz = foreCenter_xyz/count;
  foreCenter_lon_lat_rad=foreCenter_lon_lat_rad/count;
  //cout<<"foreCenter_xyz"<<foreCenter_xyz<<endl;
  cout<<"foreCenter_lon_lat_rad"<<foreCenter_lon_lat_rad<<endl;

  //TO DO: determine the bounding box on the orbital image
  Vector2 foreCenter_lon_lat;
  foreCenter_lon_lat(0) = foreCenter_lon_lat_rad(0);
  foreCenter_lon_lat(1) = foreCenter_lon_lat_rad(1);

  Vector2 backPix = backGeo.lonlat_to_pixel(foreCenter_lon_lat);
  cout<<"backCenterPix="<<backPix<<endl;
  cout <<"upsampleRatioBackImg="<<upsampleRatioBackImg<<endl;

  float leftTop_x = backPix(0) - (float)assembledWidth/(float)(2.0*upsampleRatioBackImg);
  float leftTop_y = backPix(1) - (float)assembledHeight/(float)(2.0*upsampleRatioBackImg);
  Vector2 leftTop_xy;
  leftTop_xy(0) = leftTop_x;   
  leftTop_xy(1) = leftTop_y;
  Vector2 leftTop_lonlat =  backGeo.pixel_to_lonlat(leftTop_xy); 

  H = backGeo.transform();
  //lon = H(0,0)*i+0*j + H(0,2)
  //lat = 0*i+H(1,1)*j + H(1,2)

  H(0,0) /=upsampleRatioBackImg;
  H(1,1) /=upsampleRatioBackImg;
  
  float radius = backGeo.datum().semi_major_axis();
  Vector2 point = backGeo.pixel_to_point(leftTop_xy);
  cout<<"point="<<point<<endl;

  H(0,2) = point(0);
  H(1,2) = point(1);

  cout<<"assembled transform:"<<endl;
  cout<<"H(0,0)"<<H(0,0)<<endl;
  cout<<"H(0,1)"<<H(0,1)<<endl;
  cout<<"H(0,2)"<<H(0,2)<<endl;
  cout<<"H(1,0)"<<H(1,0)<<endl;
  cout<<"H(1,1)"<<H(1,1)<<endl;
  cout<<"H(1,2)"<<H(1,2)<<endl;

  GeoReference assembledGeo = backGeo; 
  assembledGeo.set_transform(H);

  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
  cout << "num_channels = "<< num_channels << "\n";
  
  //copy the background image to the assembled image.

  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){
      
      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;

      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);
      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      
      float b_x = backPix(0);
      float b_y = backPix(1);
      
      /*
      //pixel location in the back image coordinates
      float b_y = leftTop_y + j/upsampleRatioBackImg;
      float b_x = leftTop_x + i/upsampleRatioBackImg;     
     
      //add the foreground here - START
      Vector2 backPix;
      backPix[0] = b_x;
      backPix[1] = b_y;
      */
      Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
      
      Vector2 fore_lonlat;
      fore_lonlat(0) = back_lonlat(0)-bestDeltaLonLat(0);
      fore_lonlat(1) = back_lonlat(1)-bestDeltaLonLat(1);
      
      Vector2 forePix;
      forePix = foreGeo.lonlat_to_pixel(fore_lonlat);
      
      Vector3 fore_lonlat_rad;
      fore_lonlat_rad[0]=fore_lonlat[0];
      fore_lonlat_rad[1]=fore_lonlat[1];
      fore_lonlat_rad[2]=foreImg.impl()(forePix[0], forePix[1]);
      
      //transform in xyz coords
      Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat_rad);
      
      //transform using new rotation and translation-translation only for now.
      //Vector3 transf_fore_xyz = rotation*(fore_xyz-foreCenter_xyz)+translation+foreCenter_xyz;   
      Vector3 transf_fore_xyz = fore_xyz + translation;
      
      //transform back in spherical coords
      Vector3 transf_fore_lon_lat_rad = foreGeo.datum().cartesian_to_geodetic(transf_fore_xyz);
        
      if ((isnan(foreImg.impl()(forePix[0], forePix[1]))!=FP_NAN) && (foreImg.impl()(forePix[0], forePix[1])) != 0)  {   
	assembledImg.impl()(i,j) = transf_fore_lon_lat_rad(2);
      }
      else{
	assembledImg.impl()(i,j) = backImg.impl()(b_x,b_y);
      }
  
    }
  }
  
  cout<<"translation"<<translation<<endl;

  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}

//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDRG(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    string assembledImgFilename, Vector3 translation, Matrix<float,3,3>rotation, 
                    Vector3 center, Vector2 bestDeltaLonLat)
{
 
  float maxUpsampleRatioBackImg = 40.0;
  int maxAssembledImgWidth = 4096;
  int maxAssembledImgHeight = 4096;

  cout<<"Background"<<endl;
  float numPixPerDegreeBackImg = ComputePixPerDegree(backGeo, orig_backImg.impl().cols(), orig_backImg.impl().rows(), 1);
  printf("numPixPerDegreeBackmg = %f\n", numPixPerDegreeBackImg);
   cout<<"Foreground"<<endl;
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

  if (assembledWidth > maxAssembledImgWidth){
      assembledWidth = maxAssembledImgWidth;
  }
  if (assembledHeight > maxAssembledImgHeight){
      assembledHeight = maxAssembledImgHeight;
  }
  printf("foreW = %d, foreH = %d, aW = %d, aH = %d\n", 
         orig_foreImg.impl().cols(), orig_foreImg.impl().rows(), assembledWidth, assembledHeight);
  

  Matrix<double> H;
  ImageView<typename ViewT2::pixel_type> assembledImg(assembledWidth, assembledHeight);
 
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
    = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
    = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
 
  //compute the foreground centroid
  Vector2 foreCenter_lon_lat;
  int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){ 
      if ( ((foreImg.impl()(i,j)[0]) != 0) || ((foreImg.impl()(i,j)[1]) != 0) || ((foreImg.impl()(i,j)[2]) != 0) ) {       
	
	//case of cartesian coordinates
	//get the coordinates
	Vector2 forePix(i,j);
	Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	
	//fore_lon_lat(0)=fore_lon_lat(0)+bestDeltaLonLat(0);
	//fore_lon_lat(1)=fore_lon_lat(1)+bestDeltaLonLat(1);
	
	foreCenter_lon_lat = foreCenter_lon_lat + fore_lon_lat;
	
	count++;
      }
    }
  }
  
  foreCenter_lon_lat = foreCenter_lon_lat/count;  
  cout<<"foreCenter_lon_lat="<<foreCenter_lon_lat<<endl;

  //determine the bounding box on the orbital image
  //this is correct
  Vector2 backCenterPix = backGeo.lonlat_to_pixel(foreCenter_lon_lat);
  cout<<"backCenterPix="<<backCenterPix<<endl;

  cout <<"upsampleRatioBackImg="<<upsampleRatioBackImg<<endl;
  Vector2 backLeftTopPix;

  backLeftTopPix(0) = backCenterPix(0) - (float)assembledWidth/(float)(2.0*upsampleRatioBackImg);
  backLeftTopPix(1) = backCenterPix(1) - (float)assembledHeight/(float)(2.0*upsampleRatioBackImg);
  
  cout<<"backLeftTopPix="<<backLeftTopPix<<", backWidth="<<backImg.cols()<<", backHeight="<<backImg.rows()<<endl;
  Vector2 backLeftTop_lonlat =  backGeo.pixel_to_lonlat(backLeftTopPix); 
  cout<<"backLeftTop_lonlat="<<backLeftTop_lonlat<<endl;

  H = backGeo.transform();
  //lon = H(0,0)*i+0*j + H(0,2)
  //lat = 0*i+H(1,1)*j + H(1,2)

  H(0,0) /=upsampleRatioBackImg;
  H(1,1) /=upsampleRatioBackImg;
  
  
  //new - START
  Vector2 point = backGeo.pixel_to_point(backLeftTopPix);
  cout<<"point="<<point<<endl;

  H(0,2) = point(0);
  H(1,2) = point(1);
  //new - END
  

  cout<<"assembled transform:"<<endl;
  cout<<"H(0,0)="<<H(0,0)<<endl;
  cout<<"H(0,1)="<<H(0,1)<<endl;
  cout<<"H(0,2)="<<H(0,2)<<endl;
  cout<<"H(1,0)="<<H(1,0)<<endl;
  cout<<"H(1,1)="<<H(1,1)<<endl;
  cout<<"H(1,2)="<<H(1,2)<<endl;

  GeoReference assembledGeo = backGeo; 
  assembledGeo.set_transform(H);

  cout<<"backLeftTopPix_x="<<backLeftTopPix(0)<<", backLeftTopPix_y="<<backLeftTopPix(1)<<", w="<<assembledWidth/(float)upsampleRatioBackImg<<" ,h="<<assembledHeight/(float)upsampleRatioBackImg<<endl;
   
  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
  cout << "num_channels = "<< num_channels << "\n";
  
  //copy the background image to the assembled image.
  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){

      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;

      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);
      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      
      float b_x = backPix(0);
      float b_y = backPix(1);

      Vector2 forePix = foreGeo.lonlat_to_pixel(assembledLonLat);
      float f_x = forePix(0);
      float f_y = forePix(1);
    
      if ((foreImg.impl()(f_x, f_y)[0]!=255) && (foreImg.impl()(f_x, f_y)[1]!=255) && (foreImg.impl()(f_x, f_y)[2]!=255) && 
	  (foreImg.impl()(f_x, f_y)[0]!=0) && (foreImg.impl()(f_x, f_y)[1]!=0) && (foreImg.impl()(f_x, f_y)[2]!=0)){
	  assembledImg.impl()(i,j)[0] = foreImg.impl()(f_x,f_y)[0];
	  assembledImg.impl()(i,j)[1] = foreImg.impl()(f_x,f_y)[1];
	  assembledImg.impl()(i,j)[2] = foreImg.impl()(f_x,f_y)[2];
      }
      else{
	 //cout<<"i="<<i<<", j="<<j<<", bx="<<b_x<<", b_y="<<b_y<<", backWidth="<<backImg.cols()<<", backHeight="<<backImg.rows()<<endl; 
	 assembledImg.impl()(i,j)[0] = backImg.impl()(b_x,b_y)[0];
	 assembledImg.impl()(i,j)[1] = backImg.impl()(b_x,b_y)[0];
	 assembledImg.impl()(i,j)[2] = backImg.impl()(b_x,b_y)[0];
      }
    }
  }
  
  cout<<"translation"<<translation<<endl;
  /*
  //InterpolationView<EdgeExtensionView<EdgeExtensionView<ImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> f_assembledImg
  //  = interpolate(edge_extend(assembledImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  
  //transform the foreground image and copy it over the background image
  for (int j = 0; j < orig_foreImg.impl().rows()-1; j++){
    for (int i = 0; i < orig_foreImg.impl().cols()-1; i++){
      
      if ((foreImg.impl()(i,j)[0]!=255) && (foreImg.impl()(i,j)[1]!=255) && (foreImg.impl()(i,j)[2]!=255) && 
	  (foreImg.impl()(i,j)[0]!=0) && (foreImg.impl()(i,j)[1]!=0) && (foreImg.impl()(i,j)[2]!=0)){

	//get the coordinates
	Vector2 forePix(i,j);
	Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);

	//new components - START
	fore_lon_lat(0)=fore_lon_lat(0)+bestDeltaLonLat(0);
	fore_lon_lat(1)=fore_lon_lat(1)+bestDeltaLonLat(1);
	
	//change into USGS coords
	Vector2 back_lon_lat;
	//back_lon_lat = fore_2_back_lonlat(fore_lon_lat);
        back_lon_lat = fore_lon_lat;   

        //determine their location on the assembled image
	Vector2 assembledPix = assembledGeo.lonlat_to_pixel(back_lon_lat);
        int x = assembledPix(0);
	int y = assembledPix(1);
        //cout<<"assembledImg: x="<<x<<", y="<<y<<"width= "<<assembledImg.cols()<<", height="<<assembledImg.rows()<<endl;
	      
	if ((x>0) && (y>0) && (x<assembledImg.impl().cols()) && (y<assembledImg.impl().rows())){ 
	  //overwrite the assembled image 
	  assembledImg.impl()(x,y)[0] = foreImg.impl()(i,j)[0];
	  assembledImg.impl()(x,y)[1] = foreImg.impl()(i,j)[1];
	  assembledImg.impl()(x,y)[2] = foreImg.impl()(i,j)[2];
	}
      }
    }
  }
  */
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}






#endif

