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
//#include "assembler.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;



void ComputeTileParams(int backImgWidth, int backImgHeight, 
		       GeoReference const &foreGeo, GeoReference const &backGeo, 
                       int tileSize, Vector4 paddingParams, float foreNoDataVal,
		       int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray)
{
  
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
  ComputePyramidTilingParams(backImgWidth, backImgHeight, tileSize, numPyrTiles, offsetPix, maxNumPyrLevels);
  
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

/*
//function to determine the tile properties of the area merged between the two DEMs
//the assembled DEM is derived from the georef of the background DEM
void ComputeTileParams(int orig_backImgWidth, int orig_backImgHeight, 
                       GeoReference const &foreGeo, GeoReference const &backGeo, 
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
  if (imageType == 1){//DRG
    tileSize  = assemblerParams.tileSizeDRG;
    paddingParams = assemblerParams.paddingParamsDRG; 
    foreNoDataVal = assemblerParams.foreNoDataValDRG;
  } 

  
  
  
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
  
  //int orig_backImgWidth = orig_backImg.impl().cols(); 
  //int orig_backImgHeight = orig_backImg.impl().rows();

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

//computes the parameters of the pyramid of the assembled image
//often the assembled image has the same size as the background image
//(the foreground image is a small fraction of the background) 
void ComputePyramidTilingParams(/*GeoReference const &backGeo,*/ 
                                int backImgWidth, int backImgHeight,  int tileSize,  int &numTiles, Vector2 &offsetPix, int  &numPyrLevels)
{

  cout<<"----------------------------------"<<endl;

  //int backWidth = backImg.impl().cols(); 
  //int backHeight = backImg.impl().rows();
  //cout<<"PYRAMID PARAMS: background image size: Width = "<<backWidth<<", Height = "<<backHeight<<endl;
  /*
  Vector2 backLeftTopPixel(0,0);
  Vector2 backLeftTopLonLat = backGeo.pixel_to_lonlat(backLeftTopPixel);

  Vector2 backRightBottomPixel(backImgWidth-1, backImgHeight-1);
  Vector2 backRightBottomLonLat = backGeo.pixel_to_lonlat(backRightBottomPixel);

  Vector2 assembledLeftTopLonLat = backLeftTopLonLat;
  Vector2 assembledRightBottomLonLat = backRightBottomLonLat;
  
  Vector2 assembledLeftTopPixel = backGeo.lonlat_to_pixel(assembledLeftTopLonLat);
  Vector2 assembledRightBottomPixel = backGeo.lonlat_to_pixel(assembledRightBottomLonLat);
  */
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
void ComputeSubPyramidTilingParams(GeoReference const &backGeo, 
                                   Vector4 lonlatBB, int tileSize, Vector2 offsetPix, float backUpsamplingFactor, 
                                   Vector2 &topLeftPix, Vector2 &bottomRightPix, Vector2 &topLeftTile, 
                                   Vector2 &numSubPyrTiles, int  &numSubPyrLevels)
{
  float minLon = lonlatBB(0);
  float maxLon = lonlatBB(1);
  float minLat = lonlatBB(2);
  float maxLat = lonlatBB(3);

  cout.precision(7);

  //cout<<"lonlatBB: "<<lonlatBB<<endl;
  cout<<"----------------------------------"<<endl;
  cout<<"SUBPYRAMID PARAMS: foreground region minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;

  
  //this can be removed - START
  Matrix<double> H;
  H = backGeo.transform();

  cout<<"SUBPYRAMID PARAMS: back transform:"<<endl;
  cout<<"H(0,0)="<<H(0,0)<<endl;
  cout<<"H(0,1)="<<H(0,1)<<endl;
  cout<<"H(0,2)="<<H(0,2)<<endl;
  cout<<"H(1,0)="<<H(1,0)<<endl;
  cout<<"H(1,1)="<<H(1,1)<<endl;
  cout<<"H(1,2)="<<H(1,2)<<endl;
  //this can be removed - END
  
  //cout<<"SUBPYRAMID PARAMS: lonlat bounding box: "<<lonlatBB<<endl;

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



