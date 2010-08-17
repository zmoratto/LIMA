// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

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
#include <vw/Photometry.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>

#include <stdio.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
#include "io.h"
#include "coregister.h"



void ShowFinalTrackPtsOnImage(vector<vector<LOLAShot> >trackPts, Vector<float, 6> d, 
                              string DRGFilename, string outFilename)
{
  DiskImageView<PixelRGB<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);
  
 
  vector<pointCloud> ptHere;
  int minX = DRG.cols();
  int maxX = 0;
  int minY = DRG.rows();
  int maxY = 0;
  int point_size = 7;

  for (int i = 0; i < trackPts.size(); i++){//for each track
    for(int j = 0; j < trackPts[i].size(); j++){ //for each shot in a track
      for(int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){ //for each pt in a shot 
 
	  pointCloud pt = trackPts[i][j].LOLAPt[k]; 
	  float lon = pt.coords[0];
	  float lat = pt.coords[1];
	  float rad = pt.coords[2];
	
	  Vector2 DEM_lonlat(lon, lat);
	  Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

	  int x = (int)DRG_pix[0];
	  int y = (int)DRG_pix[1];

          //compute the transformed pts
          int xd = (int)floor(d[0]*x + d[1]*y + d[2]);
          int yd = (int)floor(d[3]*x + d[4]*y + d[5]);

          if (x < minX){
	    minX = x;
          }
	  if (y < minY){
	    minY = y;
	  }
	  if (x > maxX){
	    maxX = x;
	  }
	  if (y > maxY){
	    maxY = y;
	  }
	  
	  if (xd < minX){
	    minX = xd;
	  }
	  if (yd < minY){
	       minY = yd;
	  }
	  if (xd > maxX){
	    maxX = xd;
           }
	  if (yd > maxY){
	    maxY = yd;
	  }

        
	}
      }
    }

   //printf("minX = %d, minY = %d, maxX = %d,  maxY = %d\n", minX, minY, maxX, maxY);
    
   //extend the bounding box
   maxX = maxX + 4*point_size;
   minX = minX - 4*point_size;
   maxY = maxY + 4*point_size;
   minY = minY - 4*point_size;

   //make sure the bounding box is within the image boundaries
   if (minX < 0) minX = 0;   
   if (minY < 0) minY = 0; 
   if (maxX > DRG.cols()-1) maxX = DRG.cols()-1;   
   if (maxY > DRG.rows()-1) maxY = DRG.rows()-1; 

   ImageView<PixelRGB<uint8> > DRG_crop = crop(DRG, int32(minX), int32(minY), maxX-minX+1, maxY-minY+1);
   for (int i = 0; i < trackPts.size(); i++){//for each track
    for(int j = 0; j < trackPts[i].size(); j++){ //for each shot in a track
      for(int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){ //for each pt in a shot 
 
        
	  pointCloud pt = trackPts[i][j].LOLAPt[k]; 
	  float lon = pt.coords[0];
	  float lat = pt.coords[1];
	  float rad = pt.coords[2];
	
	  Vector2 DEM_lonlat(lon, lat);
	  Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

	  int x = (int)DRG_pix[0];
	  int y = (int)DRG_pix[1];
          
    
          //if ( (x-minX > 0) && (x-minX < DRG_crop.cols()) && (y-minY > 0) && (y-minY < DRG_crop.rows())){
          //    DRG_crop(x-minX, y-minY) =  PixelRGB<uint8>(255, 0, 0);
          //}
          fill(crop(DRG_crop, int32(x) - point_size-minX, 
                              int32(y) - point_size-minY, 
                              point_size, point_size), PixelRGB<uint8>(255, 0, 0));

          //compute the transformed pts
          int xd = (int)floor(d[0]*x + d[1]*y + d[2]);
          int yd = (int)floor(d[3]*x + d[4]*y + d[5]);

        
	  //if ( (xd-minX > 0) && (xd-minX < DRG_crop.cols()) && (yd-minY > 0) && (yd-minY < DRG_crop.rows())){
          //   DRG_crop(xd-minX, yd-minY) =  PixelRGB<uint8>(0, 255, 255);
	  //}	
          fill(crop(DRG_crop, int32(xd) - point_size-minX, 
                              int32(yd) - point_size-minY, 
                              point_size, point_size), PixelRGB<uint8>(0, 255, 255));  
       
        
	}
      }
    }

  //write output image
 printf("ready to write\n");
  write_georeferenced_image(outFilename,
                            //crop(DRG, int32(minX), int32(minY), maxX-minX, maxY-minY),
                            DRG_crop,
                            DRGGeo, TerminalProgressCallback("Core","Processing:"));
}

//displays one track over the image
void ShowTrackPtsOnImage(vector<LOLAShot> trackPts, string DRGFilename, string outFilename)
{ 
  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  ImageView<PixelGray<uint8> > OutImage(DRG.cols(), DRG.rows());


  
  vector<pointCloud> ptHere;


  for (int i = 0; i < trackPts.size(); i++){//for each shot in a track
      ptHere = trackPts[i].LOLAPt;
      pointCloud centerPt  = GetPointFromIndex( ptHere, 3);
      pointCloud topPt     = GetPointFromIndex( ptHere, 2);
      pointCloud leftPt    = GetPointFromIndex( ptHere, 1);

      if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){
      
	float lon = centerPt.coords[0];
	float lat = centerPt.coords[1];
	float rad = centerPt.coords[2];

	Vector2 DEM_lonlat(lon, lat);
	Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

	int x = (int)DRG_pix[0];
	int y = (int)DRG_pix[1];

	OutImage(x,y) = 255;
      }
  }
  
  //write output image
  write_georeferenced_image(outFilename,
                            OutImage,
                            DRGGeo, TerminalProgressCallback("Core","Processing:"));
  
}


