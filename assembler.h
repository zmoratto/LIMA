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

typedef struct tilingParams
{
  float backUpsampleFactor;
  float foreUpsampleFactor;
  int back_xl;
  int back_yt;
  int back_xr;
  int back_yb;
  int horTileIndex;
  int verTileIndex;
  string filename;
  string pcFilename;
};


template <class ViewT1, class ViewT2 >
void ComputeBoundariesDEM(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
			  ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo, 
			  Vector2 bestDeltaLonLat, int tileSize, Vector2 offsetPix, Vector4 paddingParams, 
                          std::vector<struct tilingParams> &tileParamsArray)
{



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
  
  //compute centroid and boundary coordinates - START
  Vector2 foreCenterLonLat;

  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;

  int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){      
      if ((isnan(foreImg.impl()(i,j))!=FP_NAN))  {          
	  //case of cartesian coordinates
	  //get the coordinates
	  Vector2 forePix(i,j);
	  Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	  fore_lon_lat(0)=fore_lon_lat(0)+bestDeltaLonLat(0);
	  fore_lon_lat(1)=fore_lon_lat(1)+bestDeltaLonLat(1);

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
    
          foreCenterLonLat = foreCenterLonLat + fore_lon_lat;
    	  count++;
	}
    }
  }
  
  foreCenterLonLat = foreCenterLonLat/count;
  //compute centroid and boundary coordinates - END

  cout<<"---------------------------------"<<endl;
  cout<<"foreNumPixelPerDegree="<<foreNumPixPerDegree<<endl;
  cout<<"backNumPixelPerDegree="<<backNumPixPerDegree<<endl;
  
  float backUpsamplingFactor = foreNumPixPerDegree/backNumPixPerDegree;
  cout<<"backUpsamplingFactor="<<backUpsamplingFactor<<endl;

  cout<<"foreCenterLonLat"<<foreCenterLonLat<<endl;
  cout<<"minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;
  
  //compute initial pixel boundaries within the HiRISE image
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

  //determine the number of tiles
  int numHorTiles = fabs(bottomRightPix(0)-topLeftPix(0))/tileSize;
  int numVerTiles = fabs(bottomRightPix(1)-topLeftPix(1))/tileSize;
  cout<<"numHorTiles = "<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
 
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
	tilingParams thisTileParams;
   
        thisTileParams.back_xl = i - leftNumPixPadding;
        thisTileParams.back_xr = i + tileSize + rightNumPixPadding;
        thisTileParams.back_yt = j - topNumPixPadding;
        thisTileParams.back_yb = j + tileSize + bottomNumPixPadding;
        thisTileParams.horTileIndex = horTileIndex;
	thisTileParams.verTileIndex = verTileIndex;
	
	stringstream ss;
	ss<<tileParamsArray[i].horTileIndex<<"_"<<tileParamsArray[i].verTileIndex;
	string assembledDEMFilename =  /*resDir+"/"+*/"assembled_"+ss.str()+"_dem.tif";
	cout<<assembledDEMFilename<<endl;
	string assembledPCFilename =  /*resDir+"/"+*/"assembled_"+ss.str()+"_pc.txt";

        thisTileParams.filename = assembledDEMFilename;
	thisTileParams.pcFilename = assembledPCFilename;
        //if (file exists)
	thisTileParams.backUpsampleFactor = pow(2, numPyramidLevels);
	thisTileParams.foreUpsampleFactor = thisTileParams.backUpsampleFactor/backUpsamplingFactor;
	//else
	//thisTileParams.backUpsampleFactor = 1;
	//thisTileParams.foreUpsampleFactor = thisTileParams.backUpsampleFactor/backUpsamplingFactor;
	//thisTileParams.back_xl = 0;
        //thisTileParams.back_xr = tileSize + leftNumPixPadding + rightNumPixPadding;
        //thisTileParams.back_yt = 0;
        //thisTileParams.back_yb = tileSize + topNumPixPadding + bottomNumPixPadding;

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

template <class ViewT1, class ViewT2 >
void ComputeBoundariesDRG(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
			  ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo, 
			  int tileSize, Vector2 offsetPix, Vector4 paddingParams,  vector<struct tilingParams> &tileParamsArray)
{



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
  
  //compute centroid and boundary coordinates - START
  Vector2 foreCenterLonLat;

  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;

  int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){ 
    if ( ((foreImg.impl()(i,j)[0]) != 0) || ((foreImg.impl()(i,j)[1]) != 0) || ((foreImg.impl()(i,j)[2]) != 0) ) { 
       
	  //case of cartesian coordinates
	  //get the coordinates
	  Vector2 forePix(i,j);
	  Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	  fore_lon_lat(0)=fore_lon_lat(0);//+bestDeltaLonLat(0);
	  fore_lon_lat(1)=fore_lon_lat(1);//+bestDeltaLonLat(1);

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
    
          foreCenterLonLat = foreCenterLonLat + fore_lon_lat;
    	  count++;
	}
    }
  }
  
  foreCenterLonLat = foreCenterLonLat/count;
  //compute centroid and boundary coordinates - END

  cout<<"---------------------------------"<<endl;
  cout<<"foreNumPixelPerDegree="<<foreNumPixPerDegree<<endl;
  cout<<"backNumPixelPerDegree="<<backNumPixPerDegree<<endl;
  
  float backUpsamplingFactor = foreNumPixPerDegree/backNumPixPerDegree;
  cout<<"backUpsamplingFactor="<<backUpsamplingFactor<<endl;

  cout<<"foreCenterLonLat"<<foreCenterLonLat<<endl;
  cout<<"minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;
  
  //compute initial pixel boundaries within the HiRISE image
  Vector2 topLeftLonLat;
  Vector2 bottomRightLonLat;
  topLeftLonLat(0)=minLon;
  topLeftLonLat(1)=maxLat;//minLat; UBER-HACK
  bottomRightLonLat(0)=maxLon;
  bottomRightLonLat(1)=minLat;//maxLat; UBER-HACK

  cout<<"topLeftLonLat="<<topLeftLonLat<<endl;
  cout<<"bottomRightLonLat="<<bottomRightLonLat<<endl;

  //cout<<"offsetPix="<<offsetPix<<endl;
  //Vector2 padTopLeftPix = backGeo.point_to_pixel(offsetPix);
  //cout<<"padTopLeftPix="<<padTopLeftPix<<endl;

  Vector2 topLeftPix = backGeo.lonlat_to_pixel(topLeftLonLat);
  Vector2 bottomRightPix = backGeo.lonlat_to_pixel(bottomRightLonLat);
  cout<<"topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine closest tilebreak x

  topLeftPix = floor((topLeftPix+offsetPix)/(float)tileSize)*tileSize - offsetPix;
  bottomRightPix = ceil((bottomRightPix+offsetPix)/(float)tileSize)*tileSize - offsetPix;
  cout<<"adjusted to tiles: topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;
  Vector2 topLeftTile = floor((topLeftPix+offsetPix)/tileSize);
  cout<<"topLeftTile="<<topLeftTile<<endl;

  //determine the number of tiles
  int numHorTiles = fabs(bottomRightPix(0)-topLeftPix(0))/tileSize;
  int numVerTiles = fabs(bottomRightPix(1)-topLeftPix(1))/tileSize;
  cout<<"numHorTiles = "<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
 
  //pad the tiles with the N pixels where 2^N is the background upsampling ratio 
  float numPyramidLevels = ceil(log2(backUpsamplingFactor))-1;
  cout<<"numPyramidLevels="<<numPyramidLevels<<endl;

  float leftNumPixPadding = paddingParams(0);//1;
  float topNumPixPadding = paddingParams(1);//1;
  float rightNumPixPadding = paddingParams(2);//2;
  float bottomNumPixPadding = paddingParams(3);//2;
 
  int horTileIndex, verTileIndex;
  horTileIndex = topLeftTile(0); 
  for (int i = topLeftPix(0); i < bottomRightPix(0); i = i + tileSize){
    verTileIndex = topLeftTile(1);
    for (int j = topLeftPix(1); j < bottomRightPix(1); j = j + tileSize){
	
        tilingParams thisTileParams;
        thisTileParams.backUpsampleFactor = pow(2, numPyramidLevels);
	thisTileParams.foreUpsampleFactor = thisTileParams.backUpsampleFactor/backUpsamplingFactor;
        thisTileParams.back_xl = i - leftNumPixPadding;
        thisTileParams.back_xr = i + tileSize + rightNumPixPadding;
        thisTileParams.back_yt = j - topNumPixPadding;
        thisTileParams.back_yb = j + tileSize + bottomNumPixPadding;
	thisTileParams.horTileIndex = horTileIndex;
	thisTileParams.verTileIndex = verTileIndex;
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

#if 0
void ComputeTileParams(Vector2 topLeftLonLat, Vector2 bottomRightLonLat, 
                       float foreNumPixPerDegree, float backNumPixPerDegree;
                       GeoReference const &backGeo, GeoReference const &foreGeo, 
		       int tileSize, vector<struct tilingParams> &tileParamsArray)
{


  /*
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
  
  //compute centroid and boundary coordinates - START
  Vector2 foreCenterLonLat;

  float minLat =  180;
  float maxLat = -180;
  float minLon =  180;
  float maxLon = -180;

  int count = 0;
  for (int j = 0; j < foreImg.impl().rows()-1; j++){
    for (int i = 0; i < foreImg.impl().cols()-1; i++){ 
    if ( ((foreImg.impl()(i,j)[0]) != 0) || ((foreImg.impl()(i,j)[1]) != 0) || ((foreImg.impl()(i,j)[2]) != 0) ) { 
       
	  //case of cartesian coordinates
	  //get the coordinates
	  Vector2 forePix(i,j);
	  Vector2 fore_lon_lat = foreGeo.pixel_to_lonlat(forePix);
	  fore_lon_lat(0)=fore_lon_lat(0);//+bestDeltaLonLat(0);
	  fore_lon_lat(1)=fore_lon_lat(1);//+bestDeltaLonLat(1);

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
    
          foreCenterLonLat = foreCenterLonLat + fore_lon_lat;
    	  count++;
	}
    }
  }
  
  foreCenterLonLat = foreCenterLonLat/count;
  //compute centroid and boundary coordinates - END
  

  cout<<"---------------------------------"<<endl;
  cout<<"foreNumPixelPerDegree="<<foreNumPixPerDegree<<endl;
  cout<<"backNumPixelPerDegree="<<backNumPixPerDegree<<endl;
  
  float backUpsamplingFactor = foreNumPixPerDegree/backNumPixPerDegree;
  cout<<"backUpsamplingFactor="<<backUpsamplingFactor<<endl;

  cout<<"foreCenterLonLat"<<foreCenterLonLat<<endl;
  cout<<"minLon="<<minLon<<", maxLon="<<maxLon<<", minLat="<<minLat<<", maxLat="<<maxLat<<endl;
  
  //compute initial pixel boundaries within the HiRISE image
  Vector2 topLeftLonLat;
  Vector2 bottomRightLonLat;
  topLeftLonLat(0)=minLon;
  topLeftLonLat(1)=minLat;
  bottomRightLonLat(0)=maxLon;
  bottomRightLonLat(1)=maxLat;
  */

  Vector2 topLeftPix = backGeo.lonlat_to_pixel(topLeftLonLat);
  Vector2 bottomRightPix = backGeo.lonlat_to_pixel(bottomRightLonLat);
  cout<<"topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine closest tilebreak x
  topLeftPix = floor(topLeftPix/tileSize)*tileSize;
  bottomRightPix = ceil(bottomRightPix/tileSize)*tileSize;
  cout<<"adjusted to tiles: topLeftPix="<<topLeftPix<<", bottomRightPix="<<bottomRightPix<<endl;

  //determine the number of tiles
  int numHorTiles = fabs(bottomRightPix(0)-topLeftPix(0))/tileSize;
  int numVerTiles = fabs(bottomRightPix(1)-topLeftPix(1))/tileSize;
  cout<<"numHorTiles = "<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
 
  //pad the tiles with the N pixels where 2^N is the background upsampling ratio 
  float numPyramidLevels = ceil(log2(backUpsamplingFactor))-1;
  cout<<"numPyramidLevels="<<numPyramidLevels<<endl;
  
  //vector<tilingParams> tileParamsArray;
  for (int i = topLeftPix(0); i < bottomRightPix(0); i = i + tileSize){
    //STRANGE!!!!!!
      for (int j = bottomRightPix(1); j < topLeftPix(1); j = j + tileSize){
	tilingParams thisTileParams;
        thisTileParams.backUpsampleFactor = pow(2, numPyramidLevels);
	thisTileParams.foreUpsampleFactor = thisTileParams.backUpsampleFactor/backUpsamplingFactor;
        thisTileParams.back_xl = i;
        thisTileParams.back_xr = i+tileSize;
        thisTileParams.back_yt = j;
        thisTileParams.back_yb = j+tileSize;
        tileParamsArray.push_back(thisTileParams);
	
        cout<<"xl="<<thisTileParams.back_xl<<endl;
        cout<<"xr="<<thisTileParams.back_xr<<endl;
        cout<<"yt="<<thisTileParams.back_yt<<endl;
        cout<<"yb="<<thisTileParams.back_yb<<endl;
        cout<<"backUpsampleFactor="<<thisTileParams.backUpsampleFactor<<endl;
        cout<<"foreUpsampleFactor="<<thisTileParams.foreUpsampleFactor<<endl; 
       
      }
  }
  
  cout<<"---------------------------------"<<endl;
}
#endif
//overlays the foreground image on the top of the background image
//if needed it downsamples and crops the background image
//or pads it with non-data values to allow for partial ovelap between background and fore image.
template <class ViewT1, class ViewT2 >
void
ComputeAssembledDEM(ImageViewBase<ViewT1> const& orig_foreImg, GeoReference const &foreGeo,
                    ImageViewBase<ViewT2> const& orig_backImg, GeoReference const &backGeo,
                    string assembledImgFilename, string assembledAccFilename, 
                    Vector3 translation, Matrix<float,3,3>rotation, 
                    Vector3 center, Vector2 bestDeltaLonLat, struct tilingParams &tileParams)
{

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
  cout<<"****center="<<center<<endl;
  
  //determine the pixel where is the rover - START
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

      if ((isnan(foreImg.impl()(forePix[0], forePix[1]))!=FP_NAN) && (foreImg.impl()(forePix[0], forePix[1])) != 0)  {   
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
  cout<<"*******minDist="<<minDist<<", maxDist="<<maxDist<<endl;
  cout<<"translation"<<translation<<endl;

  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));

  
  //save assembled DEM accuracy as well 
  write_georeferenced_image(assembledAccFilename,
                            assembledAcc,
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
                    Vector3 center, Vector2 bestDeltaLonLat, struct tilingParams &tileParams)
{
 
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
    
      if ((foreImg.impl()(f_x, f_y)[0]!=255) && (foreImg.impl()(f_x, f_y)[1]!=255) && (foreImg.impl()(f_x, f_y)[2]!=255) && 
	  (foreImg.impl()(f_x, f_y)[0]!=0) && (foreImg.impl()(f_x, f_y)[1]!=0) && (foreImg.impl()(f_x, f_y)[2]!=0)){
	  assembledImg.impl()(i,j)[0] = foreImg.impl()(f_x,f_y)[0];
	  //assembledImg.impl()(i,j)[1] = foreImg.impl()(f_x,f_y)[1];
	  //assembledImg.impl()(i,j)[2] = foreImg.impl()(f_x,f_y)[2];
      }
      
      else{
	//cout<<"i="<<i<<", j="<<j<<", bx="<<b_x<<", b_y="<<b_y<<", backWidth="<<backImg.cols()<<", backHeight="<<backImg.rows()<<endl; 
	if((b_x<backImg.cols()) && (b_x>=0) && (b_y<backImg.rows()) && (b_y>=0)){
	  assembledImg.impl()(i,j)[0] = backImg.impl()(b_x,b_y)[0];
	  //assembledImg.impl()(i,j)[1] = backImg.impl()(b_x,b_y)[0];
	  //assembledImg.impl()(i,j)[2] = backImg.impl()(b_x,b_y)[0];
	}
      }
    }
  }
  
  //cout<<"translation"<<translation<<endl;
 
  write_georeferenced_image(assembledImgFilename,
                            assembledImg,
                            assembledGeo, TerminalProgressCallback("{Core}","Processing:"));
 
}






#endif

