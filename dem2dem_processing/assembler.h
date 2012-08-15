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
#include "pyr_tiling.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;


typedef struct AssemblerParams
{
  Vector2 deltaLonLat; //0.0001, 0.0001
  int useForeLonLatRadOffset;
  Vector2 foreLonLatOrigin;
  float foreMaxPPD;// = 2000000; 
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
  int   matchingMode;
  int   weightingMode;
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
  float deltaRad;
}; 
/*
void ComputeTileParams(int orig_backImgWidth, int orig_backImgHeight, 
		       GeoReference const &foreGeo, GeoReference const &backGeo, 
		       int tileSize, Vector4 paddingParams, float foreNoDataVal,
		       int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray);
*/
Vector4 ComputeLonLatBox(GeoReference const &geo, float noDataVal);
Vector4 ComputeLonLatBox(GeoReference const &foreGeo, float foreNoDataVal, GeoReference const &backGeo, float backNoDataVal);

void SaveAssembledPC(string DEMFilename, string assembledPCFilename);
float ComputeDEMAccuracy(GeoReference foreGeo, Vector2 forePix, float backAccuracy);

template <class ViewT1>
float ComputeTerrainResolutionPPD(ImageViewBase<ViewT1> const& img, GeoReference const &geo)
{
  int imgWidth = img.impl().cols(); 
  int imgHeight = img.impl().rows();
 
  //cout<<"imageWidth="<<imgWidth<<", imageHeight="<<imgHeight<<endl;
  
  Vector2 leftTopPixel(0,0);
  Vector2 leftTopLonLat = geo.pixel_to_lonlat(leftTopPixel);

  Vector2 rightBottomPixel(imgWidth-1, imgHeight-1);
  Vector2 rightBottomLonLat = geo.pixel_to_lonlat(rightBottomPixel);
  
  float numPixPerDegree = imgHeight/fabs(rightBottomLonLat(1)-leftTopLonLat(1));
  cout<<"TERRAIN RESOLUTION: "<<numPixPerDegree<<"PPD"<<endl;

  return numPixPerDegree;
}

//returns the average difference between two DEMs
template <class ViewT1>
float ComputeAverageDiff(ImageViewBase<ViewT1> const& foreImg, GeoReference const &foreGeo, 
			 ImageViewBase<ViewT1> const& backImg, GeoReference const &backGeo,
			 float foreNoDataVal, Vector2 offsetLonLat, struct RegistrationParams registrationParams)
{

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> interpBackImg
    = interpolate(edge_extend(backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  double delta = 0.0;
  int numValidPix = 0;
  float minDelta =  10000000;
  float maxDelta = -10000000;

  for (int i = 0; i < foreImg.impl().cols(); i=i+256){
    for (int j = 0; j < foreImg.impl().rows(); j=j+256){
      if ((isnan(foreImg.impl()(i,j))!=FP_NAN) && (foreImg.impl()(i,j)!=foreNoDataVal)){     
	Vector2 forePixel(i,j);
	Vector2 foreLonLat = foreGeo.pixel_to_lonlat(forePixel) + offsetLonLat;
        //cout<<"foreLonLat = "<<foreLonLat<<endl;
	Vector2 backPixel = backGeo.lonlat_to_pixel(foreLonLat);
	 if ((isnan(interpBackImg.impl()(backPixel(0), backPixel(1)))!=FP_NAN) && (interpBackImg.impl()(backPixel(0), backPixel(1))!=foreNoDataVal)){     
	   float thisDelta = interpBackImg.impl()(backPixel(0), backPixel(1)) - foreImg.impl()(i,j);
         
           delta = delta + thisDelta;
           if (thisDelta<minDelta){
	     minDelta = thisDelta;
           }
           if (thisDelta>maxDelta){
	     maxDelta = thisDelta;
	   }
           numValidPix++;
	 }
      } 
    }
  }  

  delta = delta/numValidPix;
  cout<<"COMPUTE AVG DIFF: delta = "<<delta<<", minDelta="<<minDelta<<", maxDelta="<<maxDelta<<endl;
  return delta;
}


//computes the  lon lat bounding box of a georeferenced image
//this function takes into account a rotation and translation of the image
//returns in foreLonLatBB the lon lat coordinates of the bounding box of the original image 
template <class ViewT1>
void ComputeLonLatBoxDEM(ImageViewBase<ViewT1> const& foreImg, GeoReference const &foreGeo, float foreNoDataVal, 
                         Vector2 deltaLonLat, struct RegistrationParams registrationParams, Vector4 &foreLonLatBB)
{

  Vector3 translation = registrationParams.translation;
  Matrix<float,3,3>rotation = registrationParams.rotation; 
  Vector3 center = registrationParams.center;
  //Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> interpForeImg
    = interpolate(edge_extend(foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

 //compute centroid and boundary coordinates of the foreground - START
  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;
 
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){      

      if ((isnan(foreImg.impl()(i,j))!=FP_NAN) && (foreImg.impl()(i,j)!=foreNoDataVal))  {     
     
	  //get the lon lat coordinates of the pixel
	  Vector2 forePix(i,j);
	  Vector2 foreLonLat = foreGeo.pixel_to_lonlat(forePix);

          //account for centroid translation
          foreLonLat(0)=foreLonLat(0) + deltaLonLat(0);
	  foreLonLat(1)=foreLonLat(1) + deltaLonLat(1);

	  //check to be within -180, 180 for longitude
	  if (foreLonLat(0) < -180){
	    foreLonLat(0) = 360 + foreLonLat(0);
	  }
	  
	  Vector3 foreLonLatRad;
          foreLonLatRad(0) = foreLonLat(0);
	  foreLonLatRad(1) = foreLonLat(1);  
	  foreLonLatRad(2) = interpForeImg.impl()(i,j);
	  
	  //change to cartesian coordinates
	  Vector3 foreXYZ = foreGeo.datum().geodetic_to_cartesian(foreLonLatRad);
	  
	  //transform using rotation and translation obtained from ICP
	  //Vector3 transfForeXYZ = rotation*(foreXYZ-center)+center+translation; 
	  Vector3 transfForeXYZ = rotation*foreXYZ+translation; 

	  //back to spherical coordinates
	  Vector3 transfForeLonLatRad = foreGeo.datum().cartesian_to_geodetic(transfForeXYZ);
  
	  //determine the location on the background image;
	  Vector2 transfForeLonLat;
	  transfForeLonLat(0) = transfForeLonLatRad(0);
	  transfForeLonLat(1) = transfForeLonLatRad(1);
	 
	  if ((isnan(transfForeLonLatRad(0))!=FP_NAN) && (isnan(transfForeLonLatRad(1))!=FP_NAN)){     
	   
	    //compute lon lat boundaries
	    if (transfForeLonLat(0) < minLon){
	      minLon = transfForeLonLat(0);
	    }
	    if (transfForeLonLat(0) > maxLon){
	      maxLon = transfForeLonLat(0);
	    }
	    if (transfForeLonLat(1) < minLat){
	      minLat = transfForeLonLat(1);
	    }
	    //cout<<transfForeLonLat(1)<<endl;
	    if (transfForeLonLat(1) > maxLat){
	      maxLat = transfForeLonLat(1);
	    }
	  }	  
      }
    }
  }
  
  foreLonLatBB(0) = minLon;
  foreLonLatBB(1) = maxLon;
  foreLonLatBB(2) = minLat;
  foreLonLatBB(3) = maxLat;
  cout<<"FOREGROUND DEM LON LAT BBox: "<<foreLonLatBB<<endl;

}


//computes the bounding box of an image
//used only when the DEM is not available
template <class ViewT1>
void ComputeLonLatBoxDRG(ImageViewBase<ViewT1> const& foreImg, GeoReference const &foreGeo, float foreNoDataVal, Vector2 deltaLonLat, Vector4 &foreLonLatBB)
{

 //compute centroid and boundary coordinates of the foreground - START
  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;
  
  //cout<<"foreImage size:"<<foreImg.impl().rows()<<", "<<foreImg.impl().rows()<<endl;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){ 
       
        if ( ((foreImg.impl()(i,j)[0]) != foreNoDataVal) || ((foreImg.impl()(i,j)[1]) != foreNoDataVal) || ((foreImg.impl()(i,j)[2]) != foreNoDataVal) ) { 

          //case of cartesian coordinates
	  Vector2 forePix(i,j);
	  Vector2 foreLonLat = foreGeo.pixel_to_lonlat(forePix);
	  foreLonLat(0)=foreLonLat(0) + deltaLonLat(0);
	  foreLonLat(1)=foreLonLat(1) + deltaLonLat(1);

	  //check to be within -180, 180 for longitude
	  if (foreLonLat(0) < -180){
	    foreLonLat(0) = 360 + foreLonLat(0);
	  }

          if (foreLonLat(0) < minLon){
	    minLon = foreLonLat(0);
          }
          if (foreLonLat(0) > maxLon){
	    maxLon = foreLonLat(0);
          }
	  if (foreLonLat(1) < minLat){
	    minLat = foreLonLat(1);
          }
          if (foreLonLat(1) > maxLat){
	    maxLat = foreLonLat(1);
          }

	}
    }
  }
  //compute boundary coordinates of the foreground- END
 
  foreLonLatBB(0) = minLon;
  foreLonLatBB(1) = maxLon;
  foreLonLatBB(2) = minLat;
  foreLonLatBB(3) = maxLat;
  cout<<"FOREGROUND DRG LON LAT BBox: "<<foreLonLatBB<<endl;

}
/*
//function to determine the tile properties of the area merged between the two DEMs
//the assembled DEM is derived from the georef of the background DEM
template <class ViewT1, class ViewT2 >
  void ComputeTileParams(ImageViewBase<ViewT1> const& orig_foreImg_t, int orig_backImgWidth, int orig_backImgHeight, GeoReference const &foreGeo,
			 ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo, 
			 //struct RegistrationParams registrationParams,
                         struct AssemblerParams assemblerParams, 
			 int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray)
{

  int tileSize;
  Vector4 paddingParams; 
  //Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;
  float foreNoDataVal;

  if (imageType == 0){//DEM
    tileSize  = assemblerParams.tileSizeDEM;
    paddingParams = assemblerParams.paddingParamsDEM; 
    foreNoDataVal = assemblerParams.foreNoDataValDEM;
  }
  if (imageType == 1){//DRG
    tileSize  = assemblerParams.tileSizeDRG;
    paddingParams = assemblerParams.paddingParamsDRG; 
    foreNoDataVal = assemblerParams.foreNoDataValDRG;
  } 

  //InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
  //  = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

  //InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
  //  = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

 
  //float foreNumPixPerDegree = ComputeTerrainResolutionPPD(orig_foreImg, foreGeo);
  //if (foreNumPixPerDegree > assemblerParams.foreMaxPPD) {foreNumPixPerDegree = assemblerParams.foreMaxPPD;}//2000000
  //float backNumPixPerDegree = ComputeTerrainResolutionPPD(orig_backImg, backGeo);

  //cout<<"COMPUTE_TILE_PARAMS: foreNumPixelPerDegree="<<foreNumPixPerDegree<<endl;
  //cout<<"COMPUTE_TILE_PARAMS: backNumPixelPerDegree="<<backNumPixPerDegree<<endl;
  //float backUpsamplingFactor = foreNumPixPerDegree/backNumPixPerDegree;
  //cout<<"COMPUTE_TILE_PARAMS: backUpsamplingFactor="<<backUpsamplingFactor<<endl;
  
  
  Matrix<double> backH;
  backH = backGeo.transform();
  cout<<"COMPUTE_TILE_PARAMS: back meter per pixel="<<backH(0,0)<<", "<<backH(1,1)<<endl;
  Matrix<double> foreH;
  foreH = foreGeo.transform();
  cout<<"COMPUTE_TILE_PARAMS: fore meter per pixel="<<foreH(0,0)<<", "<<foreH(1,1)<<endl;  
  if (foreH(0,0) < 0.015625){foreH(0,0)=0.015625;}//1/64
  float backUpsamplingFactor = fabs(backH(0,0)/foreH(0,0));
  cout<<"COMPUTE_TILE_PARAMS: backUpsamplingFactor="<<backUpsamplingFactor<<endl;
  

  int numPyrTiles;
  Vector2 offsetPix;
  int maxNumPyrLevels;
  
  int tmp_backImgWidth = orig_backImg.impl().cols(); 
  int tmp_backImgHeight = orig_backImg.impl().rows();

  ComputePyramidTilingParams(orig_backImgWidth, orig_backImgHeight, tileSize, numPyrTiles, offsetPix, maxNumPyrLevels);
  int numSubPyrLevels;
  Vector2 numSubPyrTiles;
  Vector2 topLeftPix, bottomRightPix; 
  Vector2 topLeftTile;

  ComputeSubPyramidTilingParams(backGeo, lonlatBB, tileSize, offsetPix, backUpsamplingFactor, 
                                topLeftPix, bottomRightPix, topLeftTile, numSubPyrTiles, numSubPyrLevels);

  //efectively fill in  the tiling structures - START
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
  
        thisTileParams.back_xl = i - leftNumPixPadding;
        thisTileParams.back_xr = i + tileSize + rightNumPixPadding;
        thisTileParams.back_yt = j - topNumPixPadding;
        thisTileParams.back_yb = j + tileSize + bottomNumPixPadding;
        thisTileParams.horTileIndex = horTileIndex;
	thisTileParams.verTileIndex = verTileIndex;
	
	stringstream ss;
	ss<<thisTileParams.horTileIndex<<"_"<<thisTileParams.verTileIndex;
        string assembledFilename;
	string assembledAccFilename;
	string assembledPCFilename;

        if (imageType == 0){//DEM
	  assembledFilename = "assembled_"+ss.str()+"_dem.tif";
	  cout<<"COMPUTE_TILE_PARAMS: TILE: tileDEMFilename="<<assembledFilename<<endl;
	  
	  assembledAccFilename = "assembled_"+ss.str()+"_acc.tif";
	  cout<<"COMPUTE_TILE_PARAMS: TILE: tileAccFilename="<<assembledAccFilename<<endl;
	  
	  assembledPCFilename = "assembled_"+ss.str()+"_pc.txt";
	  cout<<"COMPUTE_TILE_PARAMS: TILE: tilePCFilename="<<assembledPCFilename<<endl;
	}

        if (imageType == 1){//DRG
	  assembledFilename = "assembled_"+ss.str()+"_drg.tif";
	  cout<<"COMPUTE_TILE_PARAMS: TILE: tileDRGFilename="<<assembledFilename<<endl;
	}
       
        thisTileParams.filename = assembledFilename;
        thisTileParams.accFilename = assembledAccFilename;
	thisTileParams.pcFilename = assembledPCFilename;
	thisTileParams.backUpsampleFactor = pow(2, (float)numSubPyrLevels);
	thisTileParams.foreUpsampleFactor = thisTileParams.backUpsampleFactor/backUpsamplingFactor;

        tileParamsArray.push_back(thisTileParams);
        
        cout<<"COMPUTE_TILE_PARAMS: TILE: horTileIndex="<<thisTileParams.horTileIndex<<", verTileIndex="<<thisTileParams.verTileIndex<<endl;     
        cout<<"COMPUTE_TILE_PARAMS: TILE: xl="<<thisTileParams.back_xl<<", xr="<<thisTileParams.back_xr<<", yt="<<thisTileParams.back_yt<<", yb="<<thisTileParams.back_yb<<endl;
        cout<<"COMPUTE_TILE_PARAMS: TILE: backUpsampleFactor="<<thisTileParams.backUpsampleFactor<<endl;
        cout<<"COMPUTE_TILE_PARAMS: TILE: foreUpsampleFactor="<<thisTileParams.foreUpsampleFactor<<endl;

	verTileIndex++;
      }
    horTileIndex++;
  }
  //efectively fill in the tiling structures - END

  cout<<"---------------------------------"<<endl;
}
*/

//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
//this version removes the DEM artifacts over the foreground region
//this version replaces the old version.
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDEM(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    float foreNoDataVal,  Vector2 foreLonLatOrigin, int weightingMode, 
                    string resDirname, struct RegistrationParams registrationParams,
                    struct TilingParams &tileParams)
{

  Vector3 translation = registrationParams.translation;
  Matrix<float,3,3>rotation = registrationParams.rotation; 
  Vector3 center = registrationParams.center;
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;

  Matrix<float,3,3>inverseRotation = transpose(registrationParams.rotation);
  
  cout<<"COMPUTE_ASSEMBLED_DEM: translation="<<translation<<endl;
  cout<<"COMPUTE_ASSEMBLED_DEM: rotation="<<endl;
  PrintMatrix(rotation);
 
  float backAccuracy = 1.0;
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
 

  //set the transform of the georeference - START
  //lon = H(0,0)*i+0*j + H(0,2)
  //lat = 0*i+H(1,1)*j + H(1,2)
  H = backGeo.transform();
  H(0,0) /=upsampleRatioBackImg;
  H(1,1) /=upsampleRatioBackImg;
  //float radius = backGeo.datum().semi_major_axis();
  Vector2 point = backGeo.pixel_to_point(leftTop_xy);
  //cout<<"point="<<point<<endl;

  H(0,2) = point(0);
  H(1,2) = point(1);

  cout.precision(12);
  cout<<"COMPUTE_ASSEMBLED_DEM: assembled transform:"<<endl;
  cout<<"H(0,0)="<<H(0,0)<<endl;
  cout<<"H(0,1)="<<H(0,1)<<endl;
  cout<<"H(0,2)="<<H(0,2)<<endl;
  cout<<"H(1,0)="<<H(1,0)<<endl;
  cout<<"H(1,1)="<<H(1,1)<<endl;
  cout<<"H(1,2)="<<H(1,2)<<endl;

  GeoReference assembledGeo = backGeo; 
  assembledGeo.set_transform(H);
  //set the transform of the georeference - END


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
  //cout<<"opening"<<assembledAccFilename<<"..."<<endl; 
  //DiskImageView<float> backAcc(assembledAccFilename);  
  DiskImageView<float> backAcc = DiskImageResource::open(assembledAccFilename);  
   
  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value; 
  //determine the pixel where the rover is - START
  Vector2 camPoint;
  camPoint(0) = 0;
  camPoint(1) = 0;
  Vector2 camPix;
  camPix = foreGeo.point_to_pixel(camPoint);
  //cout<<"camPix="<<camPix<<endl;
  
  Vector2 camLonLat;
  camLonLat = foreGeo.point_to_lonlat(camPoint);
  //cout<<"camLonLat="<<camLonLat<<endl;
  //determine the pixel where is the rover - END


  //fill in the assembledImage with background values
  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){

      //cout<<"i="<<i<<", j="<<j<<", width="<<assembledWidth<<", height="<<assembledHeight<<endl;

      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;

      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);

      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      Vector2 backLonLat = backGeo.pixel_to_lonlat(backPix);

      //cout<<"test1"<<endl;

      assembledImg.impl()(i,j) = backImg.impl()(backPix[0], backPix[1]);
      //cout<<"test1.1"<<endl;
      assembledAcc(i,j) = backAcc(i,j);
      //cout<<"test1.2"<<endl;  
      
      Vector2 foreLonLat;

      float b_x = backPix(0);
      float b_y = backPix(1);
      float b_z = backImg.impl()(b_x, b_y);

      //cout<<"test1.3"<<endl;
      //determine the sub-pixel location (f_x, f_y) in the foreground image - START 
      Vector3 backLLR(assembledLonLat(0), assembledLonLat(1), b_z);
      Vector3 backXYZ = backGeo.datum().geodetic_to_cartesian(backLLR);
      Vector3 foreXYZ = inverseRotation*(backXYZ-translation);
      Vector3 foreLLR = foreGeo.datum().cartesian_to_geodetic(foreXYZ);
      //cout<<"test1.4"<<endl;

      foreLonLat(0) = foreLLR(0) - bestDeltaLonLat(0);
      foreLonLat(1) = foreLLR(1) - bestDeltaLonLat(1);
      //cout<<"test2"<<endl;

      //make sure the lon lat is within 0-180
      /*
      if (foreLonLat(0) > 180){
	//cout<<"positive: foreLonLat="<<foreLonLat<<endl;
         foreLonLat(0)=360 - foreLonLat(0);
      }
      if (foreLonLat(0) < -180){
        //cout<<"negative: foreLonLat="<<foreLonLat<<endl;
	foreLonLat(0)=360 + foreLonLat(0);
      }
      */

      //cout<<"foreLonLat="<<foreLonLat<<endl;

      Vector2 forePix = foreGeo.lonlat_to_pixel(foreLonLat);
      float f_x = forePix(0);
      float f_y = forePix(1);
      //determine the sub-pixel location (f_x, f_y) in the foreground image - END
      //cout<<"test3"<<endl;

      //if the foreground is valid, transform it using the registration params and copy it over the background.
      if ((isnan(foreImg.impl()(forePix[0], forePix[1]))!=FP_NAN) && (foreImg.impl()(forePix[0], forePix[1]) != foreNoDataVal) &&
          (forePix[0] < foreImg.cols()) && (forePix[0] >= 0) && (forePix[1] < foreImg.rows()) && (forePix[1] >= 0)){ 
	
          foreLLR(2) = foreImg.impl()(f_x,f_y);//created the interpolated foreground value
          foreXYZ = foreGeo.datum().geodetic_to_cartesian(foreLLR);	
          Vector3 transfForeXYZ = rotation*foreXYZ+translation;
          Vector3 transfForeLLR = foreGeo.datum().cartesian_to_geodetic(transfForeXYZ); 
	  
	  //cout<<"test3.1"<<endl;
          //cout<<"registrationParams.deltaRad="<<registrationParams.deltaRad<<endl;
          
          if (weightingMode == 0){
            assembledImg.impl()(i,j) = transfForeLLR(2) + registrationParams.deltaRad;
	  }
          //compute the weighted assembled value
	  if (weightingMode!=0){
	    //cout<<"test3.2"<<endl;
	      float backAccuracy = 1;//backAcc(k,l);
	      float foreAccuracy = ComputeDEMAccuracy(foreGeo, forePix, backAccuracy);
	      float foreWeight = backAccuracy/(foreAccuracy+backAccuracy);
	      assembledImg.impl()(i,j) = foreWeight*(transfForeLLR(2) + registrationParams.deltaRad) + (1-foreWeight)*backImg.impl()(b_x, b_y);
	      assembledAcc(i,j) = 1/(foreWeight*(1/foreAccuracy) + (1-foreWeight)*(1/backAccuracy)); 
	  }

      }
      
      else{ //otherwise copy the background
	if((b_x<backImg.cols()) && (b_x>=0) && (b_y<backImg.rows()) && (b_y>=0)){
	  //cout<<"test4"<<endl;
	  assembledImg.impl()(i,j) = backImg.impl()(b_x,b_y);
	}
      }
      
    }
  }

  // cout<<"test5"<<endl;
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

//NEW STYLE - WORKS WELL IN ANTARES
//there is a problem left with the DEM/DRG misalignment

//uses the DEM registartion obtained from DEMs
//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2, class ViewT3>
void
ComputeAssembledDRG(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
		    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
		    ImageViewBase<ViewT3> const& orig_backDEM, GeoReference const &backDEMGeo,
                    GeoReference &assembledGeo,
		    float foreNoDataVal, string resDirname, 
		    struct RegistrationParams registrationParams, struct TilingParams &tileParams)
{
 
  Vector3 translation = registrationParams.translation;
  Matrix<float,3,3>rotation = registrationParams.rotation; 
  Vector3 center = registrationParams.center;
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;
  string assembledImgFilename = resDirname + string("/")+tileParams.filename;
  string assembledAccFilename = resDirname + string("/")+tileParams.accFilename; 

  Matrix<float,3,3>inverseRotation = transpose(registrationParams.rotation);
   
  //new START
  int assembledWidth = (tileParams.back_xr-tileParams.back_xl)*tileParams.backUpsampleFactor;
  int assembledHeight = (tileParams.back_yb-tileParams.back_yt)*tileParams.backUpsampleFactor;
  Matrix<double> H;
  ImageView<typename ViewT2::pixel_type> assembledImg(assembledWidth, assembledHeight);

  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreImg
    = interpolate(edge_extend(orig_foreImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT2::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backImg
    = interpolate(edge_extend(orig_backImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT3::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> backDEM
    = interpolate(edge_extend(orig_backDEM.impl(),ConstantEdgeExtension()), BilinearInterpolation());
  
  Vector2 backLeftTopPix;
  backLeftTopPix(0) = tileParams.back_xl;
  backLeftTopPix(1) = tileParams.back_yt;
  
  int upsampleRatioBackImg = tileParams.backUpsampleFactor;
 
  //set the assembled georeference transform - START
  //lon = H(0,0)*i+0*j + H(0,2)
  //lat = 0*i+H(1,1)*j + H(1,2)
 
  //H = backGeo.transform();
  H = assembledGeo.transform();
  
  H(0,0) /=upsampleRatioBackImg;
  H(1,1) /=upsampleRatioBackImg;

  cout<<"COMPUTE_ASSEMBLED_DRG: backLeftTopPix="<<backLeftTopPix<<endl;
  Vector2 point = /*backGeo*/assembledGeo.pixel_to_point(backLeftTopPix);
  H(0,2) = point(0);
  H(1,2) = point(1);

  cout.precision(12);
  cout<<"COMPUTE_ASSEMBLED_DRG: assembled transform:"<<endl;
  cout<<"H(0,0)="<<H(0,0)<<endl;
  cout<<"H(0,1)="<<H(0,1)<<endl;
  cout<<"H(0,2)="<<H(0,2)<<endl;
  cout<<"H(1,0)="<<H(1,0)<<endl;
  cout<<"H(1,1)="<<H(1,1)<<endl;
  cout<<"H(1,2)="<<H(1,2)<<endl;

  //GeoReference assembledGeo = backGeo; 
  assembledGeo.set_transform(H);
  //set the assembled georeference transform - END
  
  cout<<"COMPUTE_ASSEMBLED_DRG: backLeftTopPix_x="<<backLeftTopPix(0)<<", backLeftTopPix_y="<<backLeftTopPix(1)<<", w="<<assembledWidth<<", h="<<assembledHeight<<endl;
  cout<<"COMPUTE_ASSEMBLED_DRG: bestDeltaLonLat="<<bestDeltaLonLat<<endl;
 
  //typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;
  //int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
  //cout << "num_channels = "<< num_channels << "\n";

  //copy the background image to the assembled image.
  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){
    
      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;
      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);
      Vector2 foreLonLat;

      Vector2 backDEMPix = backDEMGeo.lonlat_to_pixel(assembledLonLat);

      Vector3 backLLR(assembledLonLat(0), assembledLonLat(1), backDEM.impl()(backDEMPix(0), backDEMPix(1)));
      Vector3 backXYZ = backGeo.datum().geodetic_to_cartesian(backLLR);
      Vector3 foreXYZ = inverseRotation*(backXYZ-translation);
      Vector3 foreLLR = foreGeo.datum().cartesian_to_geodetic(foreXYZ);
      foreLonLat(0) = foreLLR(0) - bestDeltaLonLat(0);
      foreLonLat(1) = foreLLR(1) - bestDeltaLonLat(1);
      
      //foreLonLat(0) = assembledLonLat(0) - bestDeltaLonLat(0);
      //foreLonLat(1) = assembledLonLat(1) - bestDeltaLonLat(1);
      //cout<<"foreLonLat="<<foreLonLat<<endl;

      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      float b_x = backPix(0);
      float b_y = backPix(1);      

      Vector2 forePix = foreGeo.lonlat_to_pixel(foreLonLat);
      float f_x = forePix(0);
      float f_y = forePix(1);
  
      //if the foreground is valid copy it.
      if ((foreImg.impl()(f_x, f_y)[0]!=foreNoDataVal) && (foreImg.impl()(f_x, f_y)[1]!=foreNoDataVal) && (foreImg.impl()(f_x, f_y)[2]!=foreNoDataVal)){
	  assembledImg.impl()(i,j)[0] = foreImg.impl()(f_x,f_y)[0];
      }
      
      else{ //otherwise copy the background
	if((b_x<backImg.cols()) && (b_x>=0) && (b_y<backImg.rows()) && (b_y>=0)){        
	  assembledImg.impl()(i,j)[0] = backImg.impl()(b_x,b_y)[0];
	}
      }

    }
  }
  
  //cout<<"translation"<<translation<<endl;
  cout<<"COMPUTE_ASSEMBLED_DRG: assembledImgFilename="<<assembledImgFilename<<endl;
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}


//overlays the foreground image on the top of the background image
//applies NO 3D registration obtained through the DEM alignment
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2>
void
ComputeAssembledDRG(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    float foreNoDataVal, string resDirname, 
                    struct RegistrationParams registrationParams, struct TilingParams &tileParams)
{
 
  //Vector3 translation = registrationParams.translation;
  //Matrix<float,3,3>rotation = registrationParams.rotation; 
  //Vector3 center = registrationParams.center;
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;
  string assembledImgFilename = resDirname + string("/")+tileParams.filename;
  string assembledAccFilename = resDirname + string("/")+tileParams.accFilename; 

  Matrix<float,3,3>inverseRotation = transpose(registrationParams.rotation);
   

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

  cout<<"COMPUTE_ASSEMBLED_DRG: backLeftTopPix_x="<<backLeftTopPix(0)<<", backLeftTopPix_y="<<backLeftTopPix(1)<<", w="<<assembledWidth<<", h="<<assembledHeight<<endl;
  cout<<"COMPUTE_ASSEMBLED_DRG: bestDeltaLonLat="<<bestDeltaLonLat<<endl;
 
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

      Vector2 foreLonLat;

      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      float b_x = backPix(0);
      float b_y = backPix(1);

      foreLonLat(0) = assembledLonLat(0) - bestDeltaLonLat(0);
      foreLonLat(1) = assembledLonLat(1) - bestDeltaLonLat(1); 
     
      //foreLonLat(0) = assembledLonLat(0) - bestDeltaLonLat(0);
      //foreLonLat(1) = assembledLonLat(1) - bestDeltaLonLat(1);
      //cout<<"foreLonLat="<<foreLonLat<<endl;

      Vector2 forePix = foreGeo.lonlat_to_pixel(foreLonLat);
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


/*
//THIS IS A SAFE COPY 
//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
//this is very simple and cannota take into account rotations
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDRG_old(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                        ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                        float foreNoDataVal, 
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
      //cout<<"assembledLonLat="<<assembledLonLat<<endl;

      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      float b_x = backPix(0);
      float b_y = backPix(1);
      //cout<<"backPix="<<backPix<<endl;

      //cout<<"bestDeltaLonLat="<<bestDeltaLonLat<<endl;

      Vector2 foreLonLat;
      foreLonLat(0) = assembledLonLat(0) - bestDeltaLonLat(0);
      foreLonLat(1) = assembledLonLat(1) - bestDeltaLonLat(1);
      //cout<<"foreLonLat="<<foreLonLat<<endl;

      Vector2 forePix = foreGeo.lonlat_to_pixel(foreLonLat);
      float f_x = forePix(0);
      float f_y = forePix(1);

      //cout<<"forePix="<<forePix<<endl;
    
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
*/


/*
//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
//takes into account rotation and translation of the DEM
template <class ViewT1, class ViewT2, class ViewT3>
void
ComputeAssembledDRG(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
		    ImageViewBase<ViewT3> const& orig_foreDEM, GeoReference const &foreDEMGeo,
                    float foreNoDataVal, string resDirname, 
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
  InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT3::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> foreDEM
    = interpolate(edge_extend(orig_foreDEM.impl(),ConstantEdgeExtension()), BilinearInterpolation());
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
  //cout<<"point="<<point<<endl;
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
  
  //copy the background image to the assembled image - START.
  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){
 
      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;
      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);
      
      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      float b_x = backPix(0);
      float b_y = backPix(1);

      if((b_x<backImg.cols()) && (b_x>=0) && (b_y<backImg.rows()) && (b_y>=0)){
	assembledImg.impl()(i,j)[0] = backImg.impl()(b_x,b_y)[0];
      }

    }
  }
  //copy the background image to the assembled image - END.
  

  
  //copy the foreground image to the assembled image - START.
  //fill in the assembled image with foreground values
  for (int j = 0; j < foreImg.rows(); j++){
    for (int i = 0; i < foreImg.cols(); i++){
     
      Vector2 forePix;
      forePix[0] = i;
      forePix[1] = j;

      if ((isnan(foreImg.impl()(forePix[0], forePix[1])[0])!=FP_NAN) && (foreImg.impl()(forePix[0], forePix[1])[0])!=foreNoDataVal)  {  
    
        Vector2 foreLonLat = foreGeo.pixel_to_lonlat(forePix);

        //moving to background image coordinates
        Vector3 foreLonLatRad;
	foreLonLatRad(0) = foreLonLat(0) + bestDeltaLonLat(0);
	foreLonLatRad(1) = foreLonLat(1) + bestDeltaLonLat(1);
	foreLonLatRad(2) = foreDEM.impl()(forePix(0), forePix(1));

        //change to cartesian ccordinates
	Vector3 foreXYZ = assembledGeo.datum().geodetic_to_cartesian(foreLonLatRad);
      
	//transform using rotation and translation obtained from ICP
	//Vector3 transfForeXYZ = rotation*(foreXYZ-center)+center+translation; 
	Vector3 transfForeXYZ = rotation*foreXYZ+translation; 

        //back to spherical coordinates
        Vector3 transfForeLonLatRad = assembledGeo.datum().cartesian_to_geodetic(transfForeXYZ);
  
       if ((isnan(transfForeLonLatRad(0))!=FP_NAN) && (isnan(transfForeLonLatRad(1))!=FP_NAN)){     
	
	  //determine the location on the background image;
	  Vector2 transfForeLonLat;
	  transfForeLonLat(0) = transfForeLonLatRad(0);
	  transfForeLonLat(1) = transfForeLonLatRad(1);
	  
	  Vector2 assembledPix = assembledGeo.lonlat_to_pixel(transfForeLonLat);
          Vector2 backPix = backGeo.lonlat_to_pixel(transfForeLonLat);

	  int k = (int)floor(assembledPix(0));
	  int l = (int)floor(assembledPix(1));

	  if ((k>=0) && (l>=0) && (k<assembledWidth) && (l<assembledHeight)){
	    assembledImg.impl()(k,l)[0] = foreImg.impl()(i,j)[0];
	  }
	}
      }               
    }
  }
   
  //copy the foreground image into the assembled image - END.
  

  //cout<<"translation"<<translation<<endl;
  cout<<"assembledImgFilename="<<assembledImgFilename<<endl;
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}
*/

/*
//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
//this version produces DEM artifacts over the foreground region
//this version will be replaced, but so far is working very well and is the main assembler function
//KEEP THIS COPY SAFE
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDEM(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    float foreNoDataVal,  Vector2 foreLonLatOrigin, int weightingMode, 
                    string resDirname, struct RegistrationParams registrationParams,
                    struct TilingParams &tileParams)
{

  Vector3 translation = registrationParams.translation;
  Matrix<float,3,3>rotation = registrationParams.rotation; 
  Vector3 center = registrationParams.center;
  Vector2 bestDeltaLonLat = registrationParams.bestDeltaLonLat;

  float backAccuracy = 1.0;

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
  //cout<<"opening"<<assembledAccFilename<<"..."<<endl; 
  //DiskImageView<float> backAcc(assembledAccFilename);  
  DiskImageView<float> backAcc = DiskImageResource::open(assembledAccFilename);  
  
  //for (int j = 0; j < assembledHeight; j++){
  //  for (int i = 0; i < assembledWidth; i++){
  //    cout<<backAcc(i,j)<<endl;
  //  }
  // }
  
  typedef typename PixelChannelType<typename ViewT2::pixel_type>::type channel_type;

  int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value; 
  //determine the pixel where the rover is - START
  Vector2 camPoint;
  camPoint(0) = 0;
  camPoint(1) = 0;
  Vector2 camPix;
  camPix = foreGeo.point_to_pixel(camPoint);
  cout<<"camPix="<<camPix<<endl;
  
  Vector2 camLonLat;
  camLonLat = foreGeo.point_to_lonlat(camPoint);
  cout<<"camLonLat="<<camLonLat<<endl;
  //determine the pixel where is the rover - END


  //fill in the assembledImage with background values
  for (int j = 0; j < assembledHeight; j++){
    for (int i = 0; i < assembledWidth; i++){

      Vector2 assembledPix;
      assembledPix(0) = i;
      assembledPix(1) = j;

      Vector2 assembledLonLat = assembledGeo.pixel_to_lonlat(assembledPix);
      Vector2 backPix = backGeo.lonlat_to_pixel(assembledLonLat);
      Vector2 backLonLat = backGeo.pixel_to_lonlat(backPix);
      assembledImg.impl()(i,j) = backImg.impl()(backPix[0], backPix[1]);
      assembledAcc(i,j) = backAcc(i,j);

    }
  }

  //fill in the assembled image with foreground values
  for (int j = 0; j < foreImg.rows(); j++){
    for (int i = 0; i < foreImg.cols(); i++){
     
      Vector2 forePix;
      forePix[0] = i;
      forePix[1] = j;

      if ((isnan(foreImg.impl()(forePix[0], forePix[1]))!=FP_NAN) && (foreImg.impl()(forePix[0], forePix[1])) != foreNoDataVal)  {  
    
        Vector2 foreLonLat = foreGeo.pixel_to_lonlat(forePix);

        //moving to background image coordinates
        Vector3 foreLonLatRad;
	foreLonLatRad(0) = foreLonLat(0) + bestDeltaLonLat(0);
	foreLonLatRad(1) = foreLonLat(1) + bestDeltaLonLat(1);
	foreLonLatRad(2) = foreImg.impl()(forePix(0), forePix(1) );

        //change to cartesian ccordinates
	Vector3 foreXYZ = assembledGeo.datum().geodetic_to_cartesian(foreLonLatRad);
      
	//transform using rotation and translation obtained from ICP
	//Vector3 transfForeXYZ = rotation*(foreXYZ-center)+center+translation; 
   	Vector3 transfForeXYZ = rotation*foreXYZ + translation; 

        //back to spherical coordinates
        Vector3 transfForeLonLatRad = assembledGeo.datum().cartesian_to_geodetic(transfForeXYZ);
        
        transfForeLonLatRad(2) = transfForeLonLatRad(2) +registrationParams.deltaRad;   

       if ((isnan(transfForeLonLatRad(0))!=FP_NAN) && (isnan(transfForeLonLatRad(1))!=FP_NAN)){     
	
	  //determine the location on the background image;
	  Vector2 transfForeLonLat;
	  transfForeLonLat(0) = transfForeLonLatRad(0);
	  transfForeLonLat(1) = transfForeLonLatRad(1);

	  Vector2 assembledPix = assembledGeo.lonlat_to_pixel(transfForeLonLat);
          Vector2 backPix = backGeo.lonlat_to_pixel(transfForeLonLat);

	  int k = (int)floor(assembledPix(0));
	  int l = (int)floor(assembledPix(1));

	  if ((k>=0) && (l>=0) && (k<assembledWidth) && (l<assembledHeight)){
	    if (weightingMode!=0){
	      float backAccuracy = 1;//backAcc(k,l);
	      float foreAccuracy = ComputeDEMAccuracy(foreGeo, forePix, backAccuracy);
	      float foreWeight = backAccuracy/(foreAccuracy+backAccuracy);
	      
	      assembledImg.impl()(k,l) = foreWeight*(transfForeLonLatRad(2)) + (1-foreWeight)*backImg.impl()(backPix(0), backPix(1));
	      assembledAcc(k,l) = 1/(foreWeight*(1/foreAccuracy) + (1-foreWeight)*(1/backAccuracy)); 
	   
	    }
	    else{
	      assembledImg.impl()(k,l) = transfForeLonLatRad(2);
	    }
	  }
	}
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
*/


#endif


