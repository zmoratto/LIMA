// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef PYR_TILING_H
#define PYR_TILING_H

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

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

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

template <class ViewT1>
void ComputePyramidTilingParams(ImageViewBase<ViewT1> const& backImg, GeoReference const &backGeo, 
                                int tileSize,  int &numTiles, Vector2 &offsetPix, int  &numPyrLevels)
{

  cout<<"----------------------------------"<<endl;

  int backWidth = backImg.impl().cols(); 
  int backHeight = backImg.impl().rows();
 
  cout<<"PYRAMID PARAMS: image backWidth="<<backWidth<<", backHeight="<<backHeight<<endl;
  
  Vector2 backLeftTopPixel(0,0);
  Vector2 backLeftTopLonLat = backGeo.pixel_to_lonlat(backLeftTopPixel);

  Vector2 backRightBottomPixel(backWidth-1, backHeight-1);
  Vector2 backRightBottomLonLat = backGeo.pixel_to_lonlat(backRightBottomPixel);

  Vector2 assembledLeftTopLonLat = backLeftTopLonLat;
  Vector2 assembledRightBottomLonLat = backRightBottomLonLat;
  
  //cout<<"assembled image: "<<"LeftTopLonLat="<<assembledLeftTopLonLat<<", RightBottomLonLat="<<assembledRightBottomLonLat<<endl;
  Vector2 assembledLeftTopPixel = backGeo.lonlat_to_pixel(assembledLeftTopLonLat);
  Vector2 assembledRightBottomPixel = backGeo.lonlat_to_pixel(assembledRightBottomLonLat);
  //cout<<"assembled image: "<<"LeftTopPixel="<<assembledLeftTopPixel<<", RightBottomPixel="<<assembledRightBottomPixel<<endl;
  
  int assembledImgWidth = assembledRightBottomPixel(0) - assembledLeftTopPixel(0); 
  int assembledImgHeight = assembledRightBottomPixel(1) - assembledLeftTopPixel(1); 
  cout<<"PYRAMID PARAMS: assembled image size: Width = "<<assembledImgWidth<<", Height = "<<assembledImgHeight<<endl;
  //compute the assembled image size - END
  
  
  //determine the padded assembled image size - START
  float heightRatio = assembledImgHeight/(float)tileSize;
  float widthRatio = assembledImgWidth/(float)tileSize;

  float sizeRatio = heightRatio;
  if (widthRatio > sizeRatio){
      sizeRatio = widthRatio;
  }  

  //cout<<"PYRAMID PARAMS: image to tile sizeratio = "<<sizeRatio<<endl;
  numPyrLevels = ceil(log2(sizeRatio));
  cout<<"PYRAMID PARAMS: numPyrLevels="<<numPyrLevels<<endl;
 
  int padImgWidth = tileSize*pow(2, (float)numPyrLevels);    
  int padImgHeight = tileSize*pow(2, (float)numPyrLevels);
  cout<<"PYRAMID PARAMS: padded assembled image size Width="<<padImgWidth<<", Height="<<padImgHeight<<endl;
  cout<<"PYRAMID PARAMS: numTiles="<<padImgWidth/tileSize<<endl;
  cout<<"PYRAMID PARAMS: tileSize="<<tileSize<<endl;
  //determine the padded assembled image size - END


  offsetPix(0) = (padImgWidth - backImg.impl().cols())/2;
  offsetPix(1) = (padImgHeight - backImg.impl().rows())/2;
  cout<<"PYRAMID PARAMS:offsetPix="<<offsetPix<<endl;
  
  cout<<"----------------------------------"<<endl;
}

template <class ViewT1>
void ComputeSubPyramidTilingParams(ImageViewBase<ViewT1> const& backImg, GeoReference const &backGeo, Vector4 lonlatBB, 
				   int tileSize, Vector2 offsetPix, float backUpsamplingFactor, 
                                   Vector2 &topLeftPix, Vector2 &bottomRightPix, Vector2 &topLeftTile, Vector2 &numSubPyrTiles, int  &numSubPyrLevels)
{
  float minLon = lonlatBB(0);
  float maxLon = lonlatBB(1);
  float minLat = lonlatBB(2);
  float maxLat = lonlatBB(3);
  //cout<<"lonlatBB: "<<lonlatBB<<endl;
  cout<<"---------------------"<<endl;
  cout<<"SUBPYRAMID PARAMS: foreground region minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;
  
  //compute initial pixel boundaries within the HiRISE image - START
  Vector2 topLeftLonLat;
  Vector2 bottomRightLonLat;
  topLeftLonLat(0)=minLon;
  topLeftLonLat(1)=maxLat;//minLat;
  bottomRightLonLat(0)=maxLon;
  bottomRightLonLat(1)=minLat;//maxLat;

  topLeftPix = backGeo.lonlat_to_pixel(topLeftLonLat);
  bottomRightPix = backGeo.lonlat_to_pixel(bottomRightLonLat);
  cout<<"SUBPYRAMID PARAMS: BB in pixels in the original image topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine closest tilebreak x
  topLeftPix = floor((topLeftPix+offsetPix)/tileSize)*tileSize - offsetPix;
  bottomRightPix = ceil((bottomRightPix+offsetPix)/tileSize)*tileSize - offsetPix;
  cout<<"SUBPYRAMID PARAMS: BB adjusted to integer tiles. topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine the top left tile index
  topLeftTile = floor((topLeftPix+offsetPix)/tileSize);  
  cout<<"SUBPYRAMID PARAMS: topLeftTile="<<topLeftTile<<endl;
  cout<<"SUBPYRAMID PARAMS: adjusted to tiles: topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;
  //compute initial pixel boundaries within the HiRISE image - END  
  
  /*
  //determine closest tilebreak x
  topLeftPix = floor((topLeftPix+offsetPix)/tileSize)*tileSize;
  bottomRightPix = ceil((bottomRightPix+offsetPix)/tileSize)*tileSize;
  cout<<"SUBPYRAMID PARAMS: BB adjusted to integer tiles in the padded image topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine the top left tile index
  topLeftTile = floor((topLeftPix)/tileSize);  
  cout<<"SUBPYRAMID PARAMS: topLeftTile="<<topLeftTile<<endl;
  cout<<"SUBPYRAMID PARAMS: adjusted to tiles: topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;
  //compute initial pixel boundaries within the HiRISE image - END  
  */

  //determine the number of tiles for the foreground region - START
  int numHorTiles = fabs(bottomRightPix(0)-topLeftPix(0))/tileSize;
  int numVerTiles = fabs(bottomRightPix(1)-topLeftPix(1))/tileSize;
  numSubPyrTiles(0) = numHorTiles;
  numSubPyrTiles(1) = numVerTiles;
  cout<<"SUBPYRAMID PARAMS: numHorTiles = "<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
  //determine the number of tiles for the foreground region - END
  
  
  //pad the tiles with the N pixels where 2^N is the background upsampling ratio 
  numSubPyrLevels = ceil(log2(backUpsamplingFactor))-1;
  cout<<"SUBPYRAMID PARAMS: numSubPyramidLevels="<<numSubPyrLevels<<endl;
  cout<<"---------------------"<<endl;
}


#endif


