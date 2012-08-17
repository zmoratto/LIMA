// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

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

void ComputeTileParams(int backImgWidth, int backImgHeight, 
		       GeoReference const &foreGeo, GeoReference const &backGeo, 
                       int tileSize, Vector4 paddingParams, float foreNoDataVal,
		       int imageType, Vector4 lonlatBB, 
                       std::vector<struct TilingParams> &tileParamsArray)
{
  
  int numPyrLevels;
  int numPyrTiles;
  Vector2 offsetPix;
  //returns:  numPyrTiles, offsetPix and numPyrLevels
  //only the offsetPix is needed later.
  ComputePyramidTilingParams(backImgWidth, backImgHeight, tileSize, 
                             numPyrTiles, offsetPix, numPyrLevels);
 
  //needs: the offsetPix from the previous function
  //returns: backUpsamplingFactor, topLeftPix, bottomRightPix, 
  //         topLeftTile, numSubPyrTiles, numSubPyrLevels
  ComputeSubPyramidTilingParams(foreGeo, backGeo, lonlatBB, tileSize, 
                                offsetPix, paddingParams, imageType, 
                                tileParamsArray);
 
  cout<<"---------------------------------"<<endl;
}


//computes the parameters of the pyramid of the assembled image
//often the assembled image has the same size as the background image
//(the foreground image is a small fraction of the background) 
void ComputePyramidTilingParams(int backImgWidth, int backImgHeight,  int tileSize,  int &numTiles, Vector2 &offsetPix, int  &numPyrLevels)
{

 
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

  offsetPix(0) = (padImgWidth - backImgWidth)/2;
  offsetPix(1) = (padImgHeight - backImgHeight)/2;
  cout<<"PYRAMID PARAMS: offsetPix="<<offsetPix<<endl;
  
  cout<<"----------------------------------"<<endl;
}


//backGeo is in fact the assembledGeo
//topLeftPix of the corresponding tile
//botomRightPix of the corresponding tile
//index of the topLeftTile in hor and ver directions
void ComputeSubPyramidTilingParams(GeoReference const &foreGeo, GeoReference const &backGeo, 
                                   Vector4 lonlatBB, int tileSize, Vector2 offsetPix,  
                                   Vector4 paddingParams,  int imageType, 
                                   std::vector<struct TilingParams> &tileParamsArray)
{
 
  Vector2 topLeftPt, bottomRightPt; 
  Vector2 topLeftPix, bottomRightPix; 
  Vector2 topLeftTile;
  Vector2 numSubPyrTiles; 
  int  numSubPyrLevels;
  float backUpsamplingFactor;
  float foreUpsamplingFactor;  
  
  cout.precision(12);

  Matrix<double> backH;
  backH = backGeo.transform();
  cout<<"SUBPYRAMID PARAMS: back meter per pixel="<<backH(0,0)<<", "<<backH(1,1)<<endl;
  Matrix<double> foreH;
  foreH = foreGeo.transform();

  cout<<"SUBPYRAMID PARAMS: fore meter per pixel="<<foreH(0,0)<<", "<<foreH(1,1)<<endl;  
  float initForeOverBackResolutionRatio = fabs(backH(0,0)/foreH(0,0));
  if (foreH(0,0) < 0.03125){
    foreH(0,0)=0.03125;//1/32
  }

  //compute the initial resolution ratio between foreground and background
  float foreOverBackResolutionRatio = fabs(backH(0,0)/foreH(0,0));
  cout<<"SUBPYRAMID PARAMS: foreOverBackResolutionRatio="<<foreOverBackResolutionRatio<<endl;

  //compute the backUpsamplingFactor and the number of levels corresponding to a power-of-two pyramid
  numSubPyrLevels = (int)ceil(log2(foreOverBackResolutionRatio));
  backUpsamplingFactor = pow(2, (float)numSubPyrLevels);
  
  //compute the foreUpsamplingFactor
  foreUpsamplingFactor = backUpsamplingFactor/initForeOverBackResolutionRatio; 
  cout<<"SUBPYRAMID PARAMS: numSubPyrLevels="<<numSubPyrLevels<<", backUpsamplingFactor="<<backUpsamplingFactor<<", foreUpsamplingFactor="<<foreUpsamplingFactor<<endl;


  float minLon = lonlatBB(0);
  float maxLon = lonlatBB(1);
  float minLat = lonlatBB(2);
  float maxLat = lonlatBB(3);

  cout<<"SUBPYRAMID PARAMS: foreground region minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;

  Vector2 topLeftLonLat;
  Vector2 bottomRightLonLat;
  topLeftLonLat(0)=minLon;
  topLeftLonLat(1)=maxLat;//minLat;
  bottomRightLonLat(0)=maxLon;
  bottomRightLonLat(1)=minLat;//maxLat;

  //determine the point corresponding to the lonlat topleft corner
  topLeftPt = backGeo.lonlat_to_point(topLeftLonLat);
  bottomRightPt = backGeo.lonlat_to_point(bottomRightLonLat);
  cout<<"SUBPYRAMID PARAMS: BB in meters in the background image. topLeftPt="<<topLeftPt<<", bottomRightPt="<<bottomRightPt<<endl;

  //determine the pixel corresponding to the lonlat topleft corner
  topLeftPix = backGeo.lonlat_to_pixel(topLeftLonLat);
  bottomRightPix = backGeo.lonlat_to_pixel(bottomRightLonLat);
  cout<<"SUBPYRAMID PARAMS: BB in pixels in the background image. topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine closest tilebreak in pixels of the background
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
  
  //efectively fill in  the tiling structures - START
  float leftNumPixPadding = paddingParams(0);
  float topNumPixPadding = paddingParams(1);
  float rightNumPixPadding = paddingParams(2);
  float bottomNumPixPadding = paddingParams(3);
  /*
  if (imageType==0){//DEM
      leftNumPixPadding = leftNumPixPadding/4.0;
      topNumPixPadding = topNumPixPadding/4.0;
      rightNumPixPadding = rightNumPixPadding/4.0;
      bottomNumPixPadding = bottomNumPixPadding/4.0;
  }
  */
  int horTileIndex, verTileIndex;
  
  tileParamsArray.clear();
  
  horTileIndex = topLeftTile(0);
  for (int i = topLeftPix(0); i < bottomRightPix(0); i = i + tileSize){
    verTileIndex = topLeftTile(1);
    for (int j = topLeftPix(1); j < bottomRightPix(1); j = j + tileSize){
	struct TilingParams thisTileParams;
  
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
	  cout<<"SUBPYRAMID_PARAMS: TILE: tileDEMFilename="<<assembledFilename<<endl;
	  
	  assembledAccFilename = "assembled_"+ss.str()+"_acc.tif";
	  cout<<"SUBPYRAMID_PARAMS: TILE: tileAccFilename="<<assembledAccFilename<<endl;
	  
	  assembledPCFilename = "assembled_"+ss.str()+"_pc.txt";
	  cout<<"SUBPYRAMID_PARAMS: TILE: tilePCFilename="<<assembledPCFilename<<endl;
	}

        if (imageType == 1){//DRG
	  assembledFilename = "assembled_"+ss.str()+"_drg.tif";
	  cout<<"SUBPYRAMID_PARAMS: TILE: tileDRGFilename="<<assembledFilename<<endl;
	}
       
        thisTileParams.filename = assembledFilename;
        thisTileParams.accFilename = assembledAccFilename;
	thisTileParams.pcFilename = assembledPCFilename;
	thisTileParams.backUpsampleFactor = backUpsamplingFactor;
	thisTileParams.foreUpsampleFactor = foreUpsamplingFactor;

        tileParamsArray.push_back(thisTileParams);
        
        cout<<"SUBPYRAMID_PARAMS: TILE: horTileIndex="<<thisTileParams.horTileIndex<<", verTileIndex="<<thisTileParams.verTileIndex<<endl;     
        cout<<"SUBPYRAMID_PARAMS: TILE: xl="<<thisTileParams.back_xl<<", xr="<<thisTileParams.back_xr<<", yt="<<thisTileParams.back_yt<<", yb="<<thisTileParams.back_yb<<endl;
        cout<<"SUBPYRAMID_PARAMS: TILE: backUpsampleFactor="<<thisTileParams.backUpsampleFactor<<endl;
        cout<<"SUBPYRAMID_PARAMS: TILE: foreUpsampleFactor="<<thisTileParams.foreUpsampleFactor<<endl;

	verTileIndex++;
      }
    horTileIndex++;
  }
  
  
}
