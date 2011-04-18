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
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

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
                              vector<int> trackIndices, string cubFilename, string outFilename)
{

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodataVal = rsrc->nodata_read();
  cout<<"nodaval:"<<nodataVal<<endl;
  DiskImageView<PixelGray<float> > cub( rsrc );
  camera::IsisCameraModel model(cubFilename);
 

  vector<pointCloud> ptHere;
  int minX = cub.cols()-1;
  int maxX = 0;
  int minY = cub.rows()-1;
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
       
          Vector3 lon_lat_rad (lon,lat,rad*1000);
          Vector3 xyz = lon_lat_radius_to_xyz(lon_lat_rad);
          Vector2 DRG_pix = model.point_to_pixel(xyz);
          float x = DRG_pix[0];
          float y = DRG_pix[1];

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
   
   //extend the bounding box
   maxX = maxX + 4*point_size;
   minX = minX - 4*point_size;
   maxY = maxY + 4*point_size;
   minY = minY - 4*point_size;

   //make sure the bounding box is within the image boundaries
   if (minX < 0) minX = 0;   
   if (minY < 0) minY = 0; 
   if (maxX > cub.cols()-1) maxX = cub.cols()-1;   
   if (maxY > cub.rows()-1) maxY = cub.rows()-1; 
  
   printf("minX = %d, minY = %d, maxX = %d,  maxY = %d\n", minX, minY, maxX, maxY);

   ImageView<PixelGray<uint8> > DRG = apply_mask(normalize(create_mask(cub,nodataVal))*255,0);
   ImageView<PixelRGB<uint8> > DRG_crop = crop(DRG, int32(minX), int32(minY), maxX-minX+1, maxY-minY+1);
   
   //#if 0
   for (int i = 0; i < trackPts.size(); i++){//for each track
    for(int j = 0; j < trackPts[i].size(); j++){ //for each shot in a track
      for(int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){ //for each pt in a shot 
       
	//if (trackPts[i][j].valid == 1){//add here if we want just the LOLA features.
	if ((trackPts[i][j].valid == 1) && (trackPts[i][j].featurePtLOLA == 1)){//add here if we want just the LOLA features.  
            int xl, yt, w, h;

	    pointCloud pt = trackPts[i][j].LOLAPt[k]; 
	    float lon = pt.coords[0];
	    float lat = pt.coords[1];
	    float rad = pt.coords[2];
	  
            Vector3 lon_lat_rad (lon,lat,rad*1000);
            Vector3 xyz = lon_lat_radius_to_xyz(lon_lat_rad);
            Vector2 DRG_pix = model.point_to_pixel(xyz);
            float x = DRG_pix[0];
            float y = DRG_pix[1];
                  
            //make sure we are not running outside the image boundaries.
            xl = int32(x) - point_size-minX;
            yt = int32(y) - point_size-minY;
            w = point_size;
            h = point_size;
          
            if ((xl>0) && (yt>0) && (xl<DRG_crop.cols()-point_size-1) && (yt<DRG_crop.rows()-point_size-1)){
                fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(255, 0, 0));
	    }
        
	    //compute the transformed pts
	    int xd = (int)floor(d[0]*x + d[1]*y + d[2]);
	    int yd = (int)floor(d[3]*x + d[4]*y + d[5]);

            //make sure the matched point is inside the image boundaries.
	    xl = int32(xd) - point_size-minX;
            yt = int32(yd) - point_size-minY;
            w = point_size;
            h = point_size;

            if ((xl>0) && (yt>0) && (xl<DRG_crop.cols()-point_size-1) && (yt<DRG_crop.rows()-point_size-1)){
               fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(0, 0, 255));

              /*  
              if (trackPts[i][j].featurePtRefl == 1.0){
		fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(0, 255, 0));
              }
	      else{
		fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(0, 0, 255));
	      }
	      */
	    }
      
	}//valid track && LOLA feature point
      }//each point in shot
   
    }//each shot
   }
   //#endif
  //write output image
  GeoReference moonref( Datum("D_MOON"), identity_matrix<3>() );
   
  write_georeferenced_image(outFilename,
                            DRG_crop,
                            moonref, TerminalProgressCallback("Core","Processing:"));
   
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



