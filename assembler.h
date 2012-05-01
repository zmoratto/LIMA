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
#include <cmath>
#include <math.h>

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


typedef struct AssemblerParams
{
  Vector2 deltaLonLat; //0.0001, 0.0001
  int tileSizeDEM; //128
  Vector4 paddingParamsDEM; //0, 0, 1, 1
  int tileSizeDRG; //512
  Vector4 paddingParamsDRG; //1, 1, 2, 2
  float foreNoDataValDEM;
  float backNoDataValDEM;
  float foreNoDataValDRG;
  float backNoDataValDRG;
  Vector2 samplingStep;
  float maxNumStarts; 
  int matchingMode;

  Vector2 matchWindowHalfSize;
  int maxNumIter;
  float minConvThresh;
};

typedef struct RegistrationParams
{
  Vector3 translation;
  Matrix<float, 3,3 > rotation;
  Vector3 center;
  Vector2 bestDeltaLonLat;
  float error;
}; 

typedef struct TilingParams
{
  int back_xl;
  int back_yt;
  int back_xr;
  int back_yb;
  int horTileIndex;
  int verTileIndex;
  float backUpsampleFactor;
  float foreUpsampleFactor;
  string filename;
  string pcFilename;
  string accFilename;
};

Vector4 ComputeLonLatBox(GeoReference const &geo, float noDataVal);
Vector4 ComputeLonLatBox(GeoReference const &foreGeo, float foreNoDataVal, GeoReference const &backGeo, float backNoDataVal);

void SaveAssembledPC(string DEMFilename, string assembledPCFilename);

//computes the bounding box of an image
template <class ViewT1>
void ComputeLonLatBoxDEM(ImageViewBase<ViewT1> const& foreImg, GeoReference const &foreGeo, float foreNoDataVal, Vector2 deltaLonLat, Vector4 &lonlatBB)
{
  //Vector4 lonlatBB;
 //compute centroid and boundary coordinates of the foreground - START
  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;
  
  //Vector2 foreCenterLonLat;
  //int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){      
      if ((isnan(foreImg.impl()(i,j))!=FP_NAN) && (foreImg.impl()(i,j)!=foreNoDataVal))  {          
	  //case of cartesian coordinates
	  //get the coordinates
	  Vector2 forePix(i,j);
	  Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	  fore_lon_lat(0)=fore_lon_lat(0) + deltaLonLat(0);
	  fore_lon_lat(1)=fore_lon_lat(1) + deltaLonLat(1);

          if (fore_lon_lat(0) < minLon){
	    minLon = fore_lon_lat(0);
          }
          if (fore_lon_lat(0) > maxLon){
	    maxLon = fore_lon_lat(0);
          }
	  if (fore_lon_lat(1) < minLat){
	    minLat = fore_lon_lat(1);
          }
          if (fore_lon_lat(1) > maxLat){
	    maxLat = fore_lon_lat(1);
          }

          //foreCenterLonLat = foreCenterLonLat + fore_lon_lat;
    	  //count++;
	}
    }
  }
  //foreCenterLonLat = foreCenterLonLat/count;
  //cout<<"foreCenterLonLat"<<foreCenterLonLat<<endl;
  //compute centroid and boundary coordinates of the foreground- END
  lonlatBB(0) = minLon;
  lonlatBB(1) = maxLon;
  lonlatBB(2) = minLat;
  lonlatBB(3) = maxLat;
  //return lonlatBB;
}

//computes the bounding box of an image
template <class ViewT1>
void ComputeLonLatBoxDRG(ImageViewBase<ViewT1> const& foreImg, GeoReference const &foreGeo, float foreNoDataVal, Vector2 deltaLonLat, Vector4 &lonlatBB)
{
  //Vector4 lonlatBB;
 //compute centroid and boundary coordinates of the foreground - START
  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;
  
  //Vector2 foreCenterLonLat;
  //int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){      
      //if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]!=foreNoDataVal))  {     
        if ( ((foreImg.impl()(i,j)[0]) != foreNoDataVal) || ((foreImg.impl()(i,j)[1]) != foreNoDataVal) || ((foreImg.impl()(i,j)[2]) != foreNoDataVal) ) { 
	  //case of cartesian coordinates
	  //get the coordinates
	  Vector2 forePix(i,j);
	  Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	  fore_lon_lat(0)=fore_lon_lat(0) + deltaLonLat(0);
	  fore_lon_lat(1)=fore_lon_lat(1) + deltaLonLat(1);

          if (fore_lon_lat(0) < minLon){
	    minLon = fore_lon_lat(0);
          }
          if (fore_lon_lat(0) > maxLon){
	    maxLon = fore_lon_lat(0);
          }
	  if (fore_lon_lat(1) < minLat){
	    minLat = fore_lon_lat(1);
          }
          if (fore_lon_lat(1) > maxLat){
	    maxLat = fore_lon_lat(1);
          }

          //foreCenterLonLat = foreCenterLonLat + fore_lon_lat;
    	  //count++;
	}
    }
  }
  //foreCenterLonLat = foreCenterLonLat/count;
  //cout<<"foreCenterLonLat"<<foreCenterLonLat<<endl;
  //compute centroid and boundary coordinates of the foreground- END
  lonlatBB(0) = minLon;
  lonlatBB(1) = maxLon;
  lonlatBB(2) = minLat;
  lonlatBB(3) = maxLat;
  //return lonlatBB;
}

//function to determine the tile properties of the area merged between the two DEMs
//the assembled DEM is derived from the georef of the background DEM
template <class ViewT1, class ViewT2 >
void ComputeBoundaries(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
		       ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo, 
                       struct RegistrationParams registrationParams, struct AssemblerParams assemblerParams, 
                       int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray)
{

  int tileSize;
  Vector4 paddingParams; 
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;
  float foreNoDataVal;

  if (imageType == 0){//DEM
    tileSize  = assemblerParams.tileSizeDEM;
    paddingParams = assemblerParams.paddingParamsDEM; 
    foreNoDataVal = assemblerParams.foreNoDataValDEM;
  }
  if (imageType == 1){//DEM
    tileSize  = assemblerParams.tileSizeDRG;
    paddingParams = assemblerParams.paddingParamsDRG; 
    foreNoDataVal = assemblerParams.foreNoDataValDRG;
  } 

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
    = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  //compute the background an foreground upsampling factor;
  
  int foreWidth = foreImg.cols(); 
  int foreHeight = foreImg.rows();
 
  cout<<"foreWidth="<<foreWidth<<", foreHeight="<<foreHeight<<endl;

  Vector2 foreLeftTopPixel(0,0);
  Vector2 foreLeftTopLonLat = foreGeo.pixel_to_lonlat(foreLeftTopPixel);

  Vector2 foreRightBottomPixel(foreWidth-1, foreHeight-1);
  Vector2 foreRightBottomLonLat = foreGeo.pixel_to_lonlat(foreRightBottomPixel);
  
  float foreNumPixPerDegree = foreHeight/fabs(foreRightBottomLonLat(1)-foreLeftTopLonLat(1));

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
    = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  
  int backWidth = backImg.cols(); 
  int backHeight = backImg.rows();
 
  cout<<"backWidth="<<backWidth<<", backHeight="<<backHeight<<endl;

  Vector2 backLeftTopPixel(0,0);
  Vector2 backLeftTopLonLat = backGeo.pixel_to_lonlat(backLeftTopPixel);

  Vector2 backRightBottomPixel(backWidth-1, backHeight-1);
  Vector2 backRightBottomLonLat = backGeo.pixel_to_lonlat(backRightBottomPixel);
  
  float backNumPixPerDegree = backHeight/fabs(backRightBottomLonLat(1)-backLeftTopLonLat(1));

  //compute the assembled image size - START
  Vector2 assembledLeftTopLonLat = backLeftTopLonLat;
  if(fabs(foreLeftTopLonLat(0)) < fabs(assembledLeftTopLonLat(0))){
     assembledLeftTopLonLat(0) = foreLeftTopLonLat(0);
  }
  if(fabs(foreLeftTopLonLat(1)) < fabs(assembledLeftTopLonLat(1))){
     assembledLeftTopLonLat(1) = foreLeftTopLonLat(1);
  }

  Vector2 assembledRightBottomLonLat = backRightBottomLonLat;
  if(fabs(foreRightBottomLonLat(0)) > fabs(assembledRightBottomLonLat(0))){
     assembledRightBottomLonLat(0) = foreRightBottomLonLat(0);
  }
  if(fabs(foreRightBottomLonLat(1)) > fabs(assembledRightBottomLonLat(1))){
     assembledRightBottomLonLat(1) = foreRightBottomLonLat(1);
  }

  Vector2 assembledLeftTopPixel = backGeo.lonlat_to_pixel(assembledLeftTopLonLat);
  Vector2 assembledRightBottomPixel = backGeo.lonlat_to_pixel(assembledRightBottomLonLat);
  cout<<"asse: "<<"LeftTopPixel="<<assembledLeftTopPixel<<", RightBottomPixel="<<assembledRightBottomPixel<<endl;
  
  int assembledImgWidth = assembledRightBottomPixel(0) - assembledLeftTopPixel(0); 
  int assembledImgHeight = assembledRightBottomPixel(1) - assembledLeftTopPixel(1); 
  cout<<"assembled: Width = "<<assembledImgWidth<<", Height = "<<assembledImgHeight<<endl;
  //compute the assembled image size - END

  //determine the padded assembled image size - START
  float heightRatio = assembledImgHeight/(float)tileSize;//assemblerParams.tileSizeDEM;
  float widthRatio = assembledImgWidth/(float)tileSize;//assemblerParams.tileSizeDEM;
  cout<<"heightRatio = "<<heightRatio<<", widthRatio="<<widthRatio<<endl;
  float sizeRatio = heightRatio;
  if (widthRatio > sizeRatio){
      sizeRatio = widthRatio;
  }  
  cout<<"sizeratio = "<<sizeRatio<<endl;
  float maxNumPyrLevels = ceil(log2(sizeRatio));
  cout<<"maxNumPyrLevels="<<maxNumPyrLevels<<endl;

  cout<<"imgWidth="<<assembledImgWidth<<", imgHeight="<<assembledImgHeight<<endl;
  int padImgWidth = tileSize*pow(2, maxNumPyrLevels);    
  int padImgHeight = tileSize*pow(2, maxNumPyrLevels);
  cout<<"padImgWidth="<<padImgWidth<<", padImgHeight="<<padImgHeight<<endl;
  //determine the padded assembled image size - END


  Vector2 offsetPix;
  offsetPix(0) = (padImgWidth - orig_backImg.impl().cols())/2;
  offsetPix(1) = (padImgHeight - orig_backImg.impl().rows())/2;
  cout<<"offfsetPix="<<offsetPix<<endl;
 
  float minLon = lonlatBB(0);
  float maxLon = lonlatBB(1);
  float minLat = lonlatBB(2);
  float maxLat = lonlatBB(3);
  cout<<"----------------lonlatBB: "<<lonlatBB<<endl;
 
  cout<<"---------------------------------"<<endl;
  cout<<"foreNumPixelPerDegree="<<foreNumPixPerDegree<<endl;
  cout<<"backNumPixelPerDegree="<<backNumPixPerDegree<<endl;
  
  float backUpsamplingFactor = foreNumPixPerDegree/backNumPixPerDegree;
  cout<<"backUpsamplingFactor="<<backUpsamplingFactor<<endl;

  cout<<"minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;
  
  //compute initial pixel boundaries within the HiRISE image - START
  Vector2 topLeftLonLat;
  Vector2 bottomRightLonLat;
  topLeftLonLat(0)=minLon;
  topLeftLonLat(1)=maxLat;//minLat;
  bottomRightLonLat(0)=maxLon;
  bottomRightLonLat(1)=minLat;//maxLat;

  Vector2 topLeftPix = backGeo.lonlat_to_pixel(topLeftLonLat);
  Vector2 bottomRightPix = backGeo.lonlat_to_pixel(bottomRightLonLat);
  cout<<"topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine closest tilebreak x
  topLeftPix = floor((topLeftPix+offsetPix)/tileSize)*tileSize - offsetPix;
  bottomRightPix = ceil((bottomRightPix+offsetPix)/tileSize)*tileSize - offsetPix;
  //determine the top left tile index
  Vector2 topLeftTile = floor((topLeftPix+offsetPix)/tileSize);
  
  cout<<"topLeftTile="<<topLeftTile<<endl;
  cout<<"adjusted to tiles: topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;
  //compute initial pixel boundaries within the HiRISE image - END  

  //determine the number of tiles for the foreground region - START
  int numHorTiles = fabs(bottomRightPix(0)-topLeftPix(0))/tileSize;
  int numVerTiles = fabs(bottomRightPix(1)-topLeftPix(1))/tileSize;
  cout<<"numHorTiles = "<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
  //determine the number of tiles for the foreground region - END

  //pad the tiles with the N pixels where 2^N is the background upsampling ratio 
  float numPyramidLevels = ceil(log2(backUpsamplingFactor))-1;
  cout<<"numPyramidLevels="<<numPyramidLevels<<endl;


  float leftNumPixPadding = paddingParams(0);
  float topNumPixPadding = paddingParams(1);
  float rightNumPixPadding = paddingParams(2);
  float bottomNumPixPadding = paddingParams(3);

  int horTileIndex, verTileIndex;

  horTileIndex = topLeftTile(0);
  for (int i = topLeftPix(0); i < bottomRightPix(0); i = i + tileSize){
    verTileIndex = topLeftTile(1);
    for (int j = topLeftPix(1); j < bottomRightPix(1); j = j + tileSize){
	TilingParams thisTileParams;
   
        //cout<<"************horTileIndex="<<horTileIndex<<", verTileIndex="<<verTileIndex<<endl;

        thisTileParams.back_xl = i - leftNumPixPadding;
        thisTileParams.back_xr = i + tileSize + rightNumPixPadding;
        thisTileParams.back_yt = j - topNumPixPadding;
        thisTileParams.back_yb = j + tileSize + bottomNumPixPadding;
        thisTileParams.horTileIndex = horTileIndex;
	thisTileParams.verTileIndex = verTileIndex;
	
	stringstream ss;
	ss<<thisTileParams.horTileIndex<<"_"<<thisTileParams.verTileIndex;
        string assembledDEMFilename;

        if (imageType == 0){//DEM
	  assembledDEMFilename = "assembled_"+ss.str()+"_dem.tif";
	}
        if (imageType == 1){//DRG
	  assembledDEMFilename = "assembled_"+ss.str()+"_drg.tif";
	}
        cout<<assembledDEMFilename<<endl;
	
        string assembledAccFilename = "assembled_"+ss.str()+"_acc.tif";
        cout<<assembledAccFilename<<endl;

        string assembledPCFilename = "assembled_"+ss.str()+"_pc.txt";
        cout<<assembledDEMFilename<<endl;
	
    
        thisTileParams.filename = assembledDEMFilename;
        thisTileParams.accFilename = assembledAccFilename;
	thisTileParams.pcFilename = assembledPCFilename;
  
	thisTileParams.backUpsampleFactor = pow(2, numPyramidLevels);
	thisTileParams.foreUpsampleFactor = thisTileParams.backUpsampleFactor/backUpsamplingFactor;

        tileParamsArray.push_back(thisTileParams);
             
        cout<<"xl="<<thisTileParams.back_xl<<endl;
        cout<<"xr="<<thisTileParams.back_xr<<endl;
        cout<<"yt="<<thisTileParams.back_yt<<endl;
        cout<<"yb="<<thisTileParams.back_yb<<endl;
        cout<<"backUpsampleFactor="<<thisTileParams.backUpsampleFactor<<endl;
        cout<<"foreUpsampleFactor="<<thisTileParams.foreUpsampleFactor<<endl;
	cout<<"horTileIndex="<<thisTileParams.horTileIndex<<endl;
        cout<<"verTileIndex="<<thisTileParams.verTileIndex<<endl;
 
	verTileIndex++;
      }
    horTileIndex++;
  }
  
  cout<<"---------------------------------"<<endl;
}




//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDEM(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    float foreNoDataVal,  string resDirname, 
                    struct RegistrationParams registrationParams,
                    struct TilingParams &tileParams)
{

  Vector3 translation = registrationParams.translation;
  Matrix<float,3,3>rotation = registrationParams.rotation; 
  Vector3 center = registrationParams.center;
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;


  //cout<<"******* Rotation matrix "<<rotation<<endl;
  //cout<<"******* Translation vector "<<translation<<endl;
  //cout<<"******* Center "<<center<<endl;
  //cout<<"******* BestDeltaLonLat="<<bestDeltaLonLat<<endl;

  string  assembledImgFilename = resDirname + string("/")+tileParams.filename;
  string  assembledAccFilename = resDirname + string("/")+tileParams.accFilename; 
  string  assembledPCFilename  = resDirname + string("/")+tileParams.pcFilename; 

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
    = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
 
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
    = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  Vector2 leftTop_xy;
  leftTop_xy(0) = tileParams.back_xl;
  leftTop_xy(1) = tileParams.back_yt;
  int assembledWidth = (tileParams.back_xr-tileParams.back_xl)*tileParams.backUpsampleFactor;
  int assembledHeight = (tileParams.back_yb-tileParams.back_yt)*tileParams.backUpsampleFactor;
  Matrix<double> H;
  ImageView<typename ViewT2::pixel_type> assembledImg(assembledWidth, assembledHeight);
  int upsampleRatioBackImg = tileParams.backUpsampleFactor;
 
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

  //create or read the accuracy file
  ImageView<float> assembledAcc (assembledWidth, assembledHeight);
  ifstream ifile(assembledAccFilename.c_str());
  if (!ifile) {
    cout<<"writting the accuracy file"<<endl;
    for (int j = 0; j < assembledHeight; j++){
       for (int i = 0; i < assembledWidth; i++){
	 assembledAcc(i,j) = 1.0;
       }
    }
   
    //save the assembled accuracy file
    write_georeferenced_image(assembledAccFilename,
                             assembledAcc,
                             assembledGeo, TerminalProgressCallback("{Core}","Processing:"));

  }

  //open background accuracy file
  cout<<"opening"<<assembledAccFilename<<"..."<<endl; 
  DiskImageView<float> backAcc(assembledAccFilename);  
  //GeoReference backAccGeo;
  //read_georeference(backAccGeo, assembledAccFilename);
  
  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
  cout << "num_channels = "<< num_channels << "\n";
  
  //copy the background image to the assembled image.
  //cout<<"****center="<<center<<endl;
  
  //determine the pixel where the rover is - START
  Vector2 camPoint;
  camPoint(0) = 0;
  camPoint(1) = 0;
  Vector2 camPix;
  camPix = foreGeo.point_to_pixel(camPoint);
  cout<<camPix<<endl;
  Vector2 camLonLat;
  camLonLat = foreGeo.point_to_lonlat(camPoint);
  cout<<camLonLat<<endl;
  //determine the pixel where is the rover - END

  float minDist = 1000000.0;
  float maxDist = 0.0;
  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){
      
      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;

      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);
      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
     
      Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
      
      Vector2 fore_lonlat;
      fore_lonlat(0) = back_lonlat(0) - bestDeltaLonLat(0);
      fore_lonlat(1) = back_lonlat(1) - bestDeltaLonLat(1);
      
      Vector2 forePix;
      forePix = foreGeo.lonlat_to_pixel(fore_lonlat);
      
      Vector2 forePoint;
      forePoint = foreGeo.pixel_to_point(forePix);
     
      Vector3 fore_lonlat_rad;
      fore_lonlat_rad[0]=fore_lonlat[0];
      fore_lonlat_rad[1]=fore_lonlat[1];
      fore_lonlat_rad[2]=foreImg.impl()(forePix[0], forePix[1]);
   
      //transform in xyz coords
      Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat_rad);
      
      //transform using new rotation and translation-translation only for now.
      Vector3 transf_fore_xyz = rotation*(fore_xyz-center)+translation+center;   
      //Vector3 transf_fore_xyz = fore_xyz + translation;
      
      //transform back in spherical coords
      Vector3 transf_fore_lon_lat_rad = foreGeo.datum().cartesian_to_geodetic(transf_fore_xyz);
      
      float backAccuracy = backAcc(i,j);//1; //1m/pixel HIRISE accuracy

      if ((isnan(foreImg.impl()(forePix[0], forePix[1]))!=FP_NAN) && (foreImg.impl()(forePix[0], forePix[1])) != foreNoDataVal)  {   
	//determine distance from camera to current pixel
        //cout<<"forePoint="<<forePoint<<endl;
        float distToCam = sqrt(forePoint(0)*forePoint(0) + forePoint(1)*forePoint(1));
        float foreWeight = 0.5;
        
        //delta_d = 43 micro, f = 43mm, b = 0.3m
        float foreAccuracy = distToCam*distToCam/300.0; ///delta_r = r*r*delta_d/(b*f)
        foreWeight = backAccuracy/(foreAccuracy+backAccuracy);

        if (distToCam < minDist){
	  minDist = distToCam;
        }
	if (distToCam > maxDist){
	  maxDist = distToCam;
        }
        //cout<<"distToCam="<<distToCam<<endl; 
        
	//assembledImg.impl()(i,j) = transf_fore_lon_lat_rad(2);
	
	assembledImg.impl()(i,j) = foreWeight*(transf_fore_lon_lat_rad(2)) + (1-foreWeight)*backImg.impl()(backPix[0], backPix[1]);
        assembledAcc(i,j) = 1/(foreWeight*(1/foreAccuracy) + (1-foreWeight)*(1/backAccuracy)); 
      }
      else{
	assembledImg.impl()(i,j) = backImg.impl()(backPix[0], backPix[1]);
        assembledAcc(i,j) = backAccuracy;
      }
  
    }
  }
  
  //save the assembled DEM tile
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));

  
  //save assembled DEM accuracy 
  write_georeferenced_image(assembledAccFilename,
                            assembledAcc,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));

  //save assembled PC
  SaveAssembledPC(assembledImgFilename, assembledPCFilename);
  
}

//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDRG(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    float foreNoDataVal, 
                    /*string assembledImgFilename,*/
                    string resDirname, 
                    struct RegistrationParams registrationParams,
                    struct TilingParams &tileParams)
{
 
  Vector3 translation = registrationParams.translation;
  Matrix<float,3,3>rotation = registrationParams.rotation; 
  Vector3 center = registrationParams.center;
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;
  string assembledImgFilename = resDirname + string("/")+tileParams.filename;
  string assembledAccFilename = resDirname + string("/")+tileParams.accFilename; 

  //new START
  int assembledWidth = (tileParams.back_xr-tileParams.back_xl)*tileParams.backUpsampleFactor;
  int assembledHeight = (tileParams.back_yb-tileParams.back_yt)*tileParams.backUpsampleFactor;
  Matrix<double> H;
  ImageView<typename ViewT2::pixel_type> assembledImg(assembledWidth, assembledHeight);
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
    = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
    = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  Vector2 backLeftTopPix;
  backLeftTopPix(0) = tileParams.back_xl;
  backLeftTopPix(1) = tileParams.back_yt;
  int upsampleRatioBackImg = tileParams.backUpsampleFactor;
  //new END

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

  cout<<"backLeftTopPix_x="<<backLeftTopPix(0)<<", backLeftTopPix_y="<<backLeftTopPix(1)<<", w="<<assembledWidth<<", h="<<assembledHeight<<endl;
   
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
    
      //if the foreground is valid copy it.
      if ((foreImg.impl()(f_x, f_y)[0]!=foreNoDataVal) && (foreImg.impl()(f_x, f_y)[1]!=foreNoDataVal) && (foreImg.impl()(f_x, f_y)[2]!=foreNoDataVal)){
	  assembledImg.impl()(i,j)[0] = foreImg.impl()(f_x,f_y)[0];
	  //assembledImg.impl()(i,j)[1] = foreImg.impl()(f_x,f_y)[1];
	  //assembledImg.impl()(i,j)[2] = foreImg.impl()(f_x,f_y)[2];
      }
      
      else{ //otherwise copy the background
	if((b_x<backImg.cols()) && (b_x>=0) && (b_y<backImg.rows()) && (b_y>=0)){
	  assembledImg.impl()(i,j)[0] = backImg.impl()(b_x,b_y)[0];
	  //assembledImg.impl()(i,j)[1] = backImg.impl()(b_x,b_y)[0];
	  //assembledImg.impl()(i,j)[2] = backImg.impl()(b_x,b_y)[0];
	}
      }
    }
  }
  
  //cout<<"translation"<<translation<<endl;
  cout<<"assembledImgFilename="<<assembledImgFilename<<endl;
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}

#endif


