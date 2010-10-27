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
#include "tracks.h"


//displays the original tracks(red) and transformed tracks(cyan) over the image
void ShowFinalTrackPtsOnImage(vector<vector<LOLAShot> >trackPts, Vector<float, 6> d, 
                              vector<int> trackIndices, string DRGFilename, string outFilename)
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

  //for (int i = 0; i < trackPts.size(); i++){//for each track
 for (int i = 0; i < trackIndices.size(); i++){//for each selected track 
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
       
	if (trackPts[i][j].valid == 1){
	    pointCloud pt = trackPts[i][j].LOLAPt[k]; 
	    float lon = pt.coords[0];
	    float lat = pt.coords[1];
	    float rad = pt.coords[2];
	
	    Vector2 DEM_lonlat(lon, lat);
	    Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);
	  
	    int x = (int)DRG_pix[0];
	    int y = (int)DRG_pix[1];
          

	    fill(crop(DRG_crop, int32(x) - point_size-minX, 
                              int32(y) - point_size-minY, 
                              point_size, point_size), PixelRGB<uint8>(255, 0, 0));

	    //compute the transformed pts
	    int xd = (int)floor(d[0]*x + d[1]*y + d[2]);
	    int yd = (int)floor(d[3]*x + d[4]*y + d[5]);

	    fill(crop(DRG_crop, int32(xd) - point_size-minX, 
		      int32(yd) - point_size-minY, 
		      point_size, point_size), PixelRGB<uint8>(0, 255, 255));  
	}
        
	}
      }
    }

  //write output image

  write_georeferenced_image(outFilename,
                            DRG_crop,
                            DRGGeo, TerminalProgressCallback("Core","Processing:"));
}



//displays the LOLA tracks in an image format on a black background
//quick space efficient way to visualize Lidar data without storing the real image but a black background. 
void MakeGrid(vector<vector<LOLAShot> >trackPts, int numVerPts, int numHorPts, string DEMFilename, vector<int> trackIndices)
{
  int l, m, n;  
  ImageView<PixelGray<float> > DEMImage(numHorPts, numVerPts);
  GeoReference DEMgeo;

  
  if (trackIndices.size() == 1){
      //string outDEMFilename;
      //char* outDEMFilename_char = new char[500];
      //sprintf (outDEMFilename_char, "../results/dem_%d.tiff", k);
      //outDEMFilename = std::string(outDEMFilename_char);
  }

  Vector4 coords = FindMinMaxLat(trackPts);

  printf("minLat=%f, maxLat=%f, minLon=%f maxLon=%f\n", coords(0), coords(1), coords(2), coords(3));

  float minLat = coords(0); //this causes seg-fault
  float maxLat = coords(1); 
  float minLon = coords(2); 
  float maxLon = coords(3);

  float lonDelta = (maxLon-minLon)/numHorPts;
  float latDelta = (maxLat-minLat)/numVerPts;

  //printf("lonDelta = %f, latDelta = %f\n", lonDelta, latDelta);
  //init the DEM
  for (l = 0; l < numHorPts; l++){
    for (m = 0; m < numVerPts; m++){
      DEMImage(l, m) = 0.0;
    }
  }
  //fill the DEM

  printf("numTracks = %d\n", (int)(trackIndices.size()));
  for (int k = 0; k < trackIndices.size();k++){
    int trackIndex = trackIndices[k];
    printf("trackIndex = %d\n", trackIndex);
    for (n = 0; n < trackPts[trackIndex].size(); n++){ 
      for (int s = 0; s < trackPts[trackIndex][n].LOLAPt.size(); s++){

        float lon_index = (trackPts[trackIndex][n].LOLAPt[s].coords[0] - minLon)/lonDelta;
        float lat_index = (trackPts[trackIndex][n].LOLAPt[s].coords[1] - minLat)/latDelta;
        l = (int)floor(lon_index);
        m = (int)floor(lat_index);

        if ((m < numVerPts) && (l<numHorPts)){ 
          DEMImage(l, m) = trackPts[trackIndex][n].LOLAPt[s].coords[2]; 
        }
        else{
          printf("Error\n");
          printf("l = %d, m = %d, numHorPts = %d, numVerPts = %d\n", l, m, numHorPts, numVerPts);
        }
      }
    }
  }
  write_georeferenced_image(DEMFilename, 
      DEMImage,
      DEMgeo, TerminalProgressCallback("{Core}","Processing:"));

}



