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


/*
void ComputeTileParams(int orig_backImgWidth, int orig_backImgHeight, 
                       GeoReference const &foreGeo, GeoReference const &backGeo, 
		       struct RegistrationParams registrationParams, struct AssemblerParams assemblerParams, 
		       int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray);
*/
void ComputePyramidTilingParams(int backImgWidth, int backImgHeight,  int tileSize,  
                                int &numTiles, Vector2 &offsetPix, int  &numPyrLevels);
void ComputeSubPyramidTilingParams(GeoReference const &backGeo, Vector4 lonlatBB, int tileSize, Vector2 offsetPix, float backUpsamplingFactor, 
                                   Vector2 &topLeftPix, Vector2 &bottomRightPix, Vector2 &topLeftTile, 
                                   Vector2 &numSubPyrTiles, int  &numSubPyrLevels);
void ComputeTileParams(int orig_backImgWidth, int orig_backImgHeight, 
		       GeoReference const &foreGeo, GeoReference const &backGeo, 
		       int tileSize, Vector4 paddingParams, float foreNoDataVal,
		       int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray);
/*
//computes the parameters of the pyramid of the assembled image
//often the assembled image has the same size as the background image
//(the foreground image is a small fraction of the background) 
template <class ViewT1>
void ComputePyramidTilingParams(ImageViewBase<ViewT1> const& backImg_t, GeoReference const &backGeo, 
                                int backImgWidth, int backImgHeight,  int tileSize,  int &numTiles, Vector2 &offsetPix, int  &numPyrLevels)
{

  cout<<"----------------------------------"<<endl;

  //int backWidth = backImg.impl().cols(); 
  //int backHeight = backImg.impl().rows();
  //cout<<"PYRAMID PARAMS: background image size: Width = "<<backWidth<<", Height = "<<backHeight<<endl;
  
  Vector2 backLeftTopPixel(0,0);
  Vector2 backLeftTopLonLat = backGeo.pixel_to_lonlat(backLeftTopPixel);

  Vector2 backRightBottomPixel(backImgWidth-1, backImgHeight-1);
  Vector2 backRightBottomLonLat = backGeo.pixel_to_lonlat(backRightBottomPixel);

  Vector2 assembledLeftTopLonLat = backLeftTopLonLat;
  Vector2 assembledRightBottomLonLat = backRightBottomLonLat;
  
  Vector2 assembledLeftTopPixel = backGeo.lonlat_to_pixel(assembledLeftTopLonLat);
  Vector2 assembledRightBottomPixel = backGeo.lonlat_to_pixel(assembledRightBottomLonLat);

  //int assembledImgWidth = backWidth;//assembledRightBottomPixel(0) - assembledLeftTopPixel(0) + 1; 
  //int assembledImgHeight = backHeight;//assembledRightBottomPixel(1) - assembledLeftTopPixel(1) + 1;

  int assembledImgWidth = backImgWidth;
  int assembledImgHeight = backImgHeight; 
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


  //offsetPix(0) = (padImgWidth - backImg.impl().cols())/2;
  //offsetPix(1) = (padImgHeight - backImg.impl().rows())/2;

  offsetPix(0) = (padImgWidth - backImgWidth)/2;
  offsetPix(1) = (padImgHeight - backImgHeight)/2;
  cout<<"PYRAMID PARAMS: offsetPix="<<offsetPix<<endl;
  
  cout<<"----------------------------------"<<endl;
}

//note that the backImg_t is not used here and can be removed
template <class ViewT1>
void ComputeSubPyramidTilingParams(ImageViewBase<ViewT1> const& backImg_t, GeoReference const &backGeo, 
                                   Vector4 lonlatBB, int tileSize, Vector2 offsetPix, float backUpsamplingFactor, 
                                   Vector2 &topLeftPix, Vector2 &bottomRightPix, Vector2 &topLeftTile, 
                                   Vector2 &numSubPyrTiles, int  &numSubPyrLevels)
{
  float minLon = lonlatBB(0);
  float maxLon = lonlatBB(1);
  float minLat = lonlatBB(2);
  float maxLat = lonlatBB(3);

  //cout<<"lonlatBB: "<<lonlatBB<<endl;
  cout<<"----------------------------------"<<endl;
  cout<<"SUBPYRAMID PARAMS: foreground region minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;

  //this can be removed - START
  Matrix<double> H;
  H = backGeo.transform();
  cout.precision(7);
  cout<<"SUBPYRAMID PARAMS: back transform:"<<endl;
  cout<<"H(0,0)="<<H(0,0)<<endl;
  cout<<"H(0,1)="<<H(0,1)<<endl;
  cout<<"H(0,2)="<<H(0,2)<<endl;
  cout<<"H(1,0)="<<H(1,0)<<endl;
  cout<<"H(1,1)="<<H(1,1)<<endl;
  cout<<"H(1,2)="<<H(1,2)<<endl;
  //this can be removed - END

  //compute initial pixel boundaries within the HiRISE image - START
  Vector2 topLeftLonLat;
  Vector2 bottomRightLonLat;
  topLeftLonLat(0)=minLon;
  topLeftLonLat(1)=maxLat;//minLat;
  bottomRightLonLat(0)=maxLon;
  bottomRightLonLat(1)=minLat;//maxLat;

  topLeftPix = backGeo.lonlat_to_pixel(topLeftLonLat);
  bottomRightPix = backGeo.lonlat_to_pixel(bottomRightLonLat);
  cout<<"SUBPYRAMID PARAMS: BB in pixels in the background image. topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine closest tilebreak x
  topLeftPix = floor((topLeftPix+offsetPix)/tileSize)*tileSize - offsetPix;
  bottomRightPix = ceil((bottomRightPix+offsetPix)/tileSize)*tileSize - offsetPix;
  cout<<"SUBPYRAMID PARAMS: BB in pixels in the background image adjusted to tiles. topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine the top left tile index
  topLeftTile = floor((topLeftPix+offsetPix)/tileSize);  
  cout<<"SUBPYRAMID PARAMS: topLeftTile="<<topLeftTile<<endl;

  //determine the number of tiles for the foreground region - START
  numSubPyrTiles(0) = fabs(bottomRightPix(0)-topLeftPix(0))/tileSize; //numHorTiles
  numSubPyrTiles(1) = fabs(bottomRightPix(1)-topLeftPix(1))/tileSize; //numVerTiles
  cout<<"SUBPYRAMID PARAMS: numHorTiles = "<<numSubPyrTiles(0)<<", numVerTiles="<<numSubPyrTiles(1)<<endl;
  //determine the number of tiles for the foreground region - END
  
  //pad the tiles with the N pixels where 2^N is the background upsampling ratio 
  numSubPyrLevels = ceil(log2(backUpsamplingFactor))-1;
  cout<<"SUBPYRAMID PARAMS: numSubPyramidLevels="<<numSubPyrLevels<<endl;
  cout<<"----------------------------------"<<endl;
 
  
  
}
*/

#endif


