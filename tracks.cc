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

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
#include "io.h"
#include "match.h"
#include "coregister.h"
#include "display.h"


Vector2 ComputeMinMaxValuesFromCub(string cubFilename)
{
  Vector2 minmax;
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodataVal = rsrc->nodata_read();
  cout<<"nodaval:"<<nodataVal<<endl;
  DiskImageView<PixelGray<float> > isis_view( rsrc );
  int width = isis_view.cols();
  int height = isis_view.rows();
  float minVal = 100000000.0;
  float maxVal = -100000000.0;
  for (int i = 0; i < height; i++){
    for (int j = 0; j < width; j++){
      
      if ((isis_view(j,i) < minVal) && (isis_view(j,i) > nodataVal)){
	minVal = isis_view(j,i);
      }
      if ((isis_view(j,i) > maxVal) && (isis_view(j,i) > nodataVal)){
	maxVal = isis_view(j,i);
      }
    }
  }
  minmax(0) = minVal;
  minmax(1) = maxVal;
  cout<<"min="<<minVal<<", max="<<maxVal<<endl;

  return minmax;
}

void 
GetAllPtsFromCub(vector<vector<LOLAShot > > &trackPts, string cubFilename)
{

  vector<pointCloud> LOLAPts;
  
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodata_value = rsrc->nodata_read();
 
  DiskImageView<PixelGray<float> > isis_view( rsrc );
  int width = isis_view.cols();
  int height = isis_view.rows();
  camera::IsisCameraModel model(cubFilename);

  //calculate the min max value of the image
  Vector2 minmax = ComputeMinMaxValuesFromCub(cubFilename);

  ImageViewRef<float> interpImg;
  interpImg = pixel_cast<float>(interpolate(edge_extend(isis_view,
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );
  
  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      
      trackPts[k][i].valid = 1; 

      LOLAPts = trackPts[k][i].LOLAPt;
      trackPts[k][i].imgPt.resize(LOLAPts.size());

      for (int j = 0; j < LOLAPts.size(); j++){
          
	    float lon = LOLAPts[j].coords[0];
	    float lat = LOLAPts[j].coords[1];
	    float rad = LOLAPts[j].coords[2];
            
            Vector3 lon_lat_rad (lon,lat,rad*1000);
            Vector3 xyz = lon_lat_radius_to_xyz(lon_lat_rad);
            Vector2 cub_pix = model.point_to_pixel(xyz);
            float x = cub_pix[0];
            float y = cub_pix[1];
            //check that (x,y) are within the image boundaries
	    if ((x>=0) && (y>=0) && (x<width) && (y<height)){//valid position  
              //check for valid data as well
              if (interpImg(x, y)>nodata_value){//valid values
		 trackPts[k][i].imgPt[j].val = interpImg(x, y) - minmax(0);
	         trackPts[k][i].imgPt[j].x = cub_pix[0];
	         trackPts[k][i].imgPt[j].y = cub_pix[1];
	      }
              else{//invalidate the point
                 trackPts[k][i].valid = 0;
              }
	    }
            else{ //invalidate the point  
                 trackPts[k][i].valid = 0; 
            }
	
      }
      //check for valid shot with 3 points to compute reflectance
      pointCloud centerPt  = GetPointFromIndex( LOLAPts, 3);
      pointCloud topPt     = GetPointFromIndex( LOLAPts, 2);
      pointCloud leftPt    = GetPointFromIndex( LOLAPts, 1);
      if ((centerPt.s == -1) || (topPt.s == -1) || (leftPt.s == -1) || (LOLAPts.size() >5)){//invalid LOLA shot
          trackPts[k][i].valid = 0; 
      }

    }//i  
  }//k
 
}

vector<float> GetTrackPtsByID(vector<LOLAShot> trackPts, int ID)
{
  vector<float> pts;
  for (int i = 0; i < trackPts.size(); i++){
    for (int k = 0; k < trackPts[i].LOLAPt.size(); k++){

      float rad = trackPts[i].LOLAPt[k].coords[2];
      int id = trackPts[i].LOLAPt[k].s;

      if (id == ID){
        pts.push_back(rad);
      }

    }
  }

  return pts;
}

//computes the scale factor for all tracks at once
float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts)
{
  float nominator = 0.0;
  int numValidPts = 0;
  float scaleFactor = 1;

  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      if ((trackPts[k][i].valid == 1) && (trackPts[k][i].reflectance != 0) &&(trackPts[k][i].reflectance != -1)){//valid track and non-zero reflectance

        //update the nominator for the center point
        for (int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){
          if (trackPts[k][i].LOLAPt[j].s == 3){ 
            nominator = nominator + trackPts[k][i].imgPt[j].val/trackPts[k][i].reflectance;
          }
        }
        //update the denominator
        numValidPts++;
      }
    }
  }

  if (numValidPts != 0){ 
    scaleFactor = nominator/numValidPts;
  }
  return scaleFactor;
}

void ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks, ModelParams modelParams,  CoregistrationParams coregistrationParams)
{

  vector<pointCloud> LOLAPts;


  GlobalParams globalParams;
  globalParams.reflectanceType = coregistrationParams.reflectanceType;
  globalParams.slopeType = 1;
  globalParams.shadowThresh = 40;
  globalParams.albedoInitType = 1;
  globalParams.exposureInitType = 1;


  for (int k = 0; k < allTracks.size(); k++ ){
    for (int i = 0; i < allTracks[k].size(); i++){
      LOLAPts = allTracks[k][i].LOLAPt;
      pointCloud centerPt = GetPointFromIndex(LOLAPts, 3);
      pointCloud topPt = GetPointFromIndex(LOLAPts, 2);
      pointCloud leftPt = GetPointFromIndex(LOLAPts, 1);

      if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1) && (allTracks[k][i].valid == 1)){
        Datum moon;
        moon.set_well_known_datum("D_MOON");

        centerPt.coords[2] = (centerPt.coords[2]-1737.4)*1000;
        topPt.coords[2] = (topPt.coords[2]-1737.4)*1000;
        leftPt.coords[2] = (leftPt.coords[2]-1737.4)*1000;

        Vector3 xyz = moon.geodetic_to_cartesian(centerPt.coords);
        Vector3 xyzTop = moon.geodetic_to_cartesian(topPt.coords);
        Vector3 xyzLeft = moon.geodetic_to_cartesian(leftPt.coords);
        Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);

        allTracks[k][i].reflectance = ComputeReflectance(normal, xyz, modelParams, globalParams);
      }
      else{
        allTracks[k][i].reflectance = -1;
      }
    }//i
  }//k

}

pointCloud GetPointFromIndex(vector<pointCloud> const &  LOLAPts, int index)
{
  pointCloud pt;
  pt.s = -1;//invalid pointCloud
  for(int i = 0;i < LOLAPts.size(); i++){
    if (LOLAPts[i].s == index){
      return LOLAPts[i];
    }
  }

  return pt;
}


//this function will be similar to GetAllPtsFromImage: 
//TO DO: rename to GetAllPtsFromDEM(vector<vector<LOLAShot > > &trackPts,  ImageViewBase<ViewT> const& DEM, GeoReference const &DEMGeo)
vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID)
{
  DiskImageView<PixelGray<float> >   DEM(DEMFilename);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, DEMFilename);

  vector<float> demPts;
  demPts.resize(trackPts.size());

  ImageViewRef<PixelGray<float>  >  interpDEM = interpolate(edge_extend(DEM.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  int index = 0;
  for (int i = 0; i < trackPts.size(); i++){
    demPts[index] = -1;
    for(int j = 0; j < trackPts[i].LOLAPt.size(); j++){

      float lon = trackPts[i].LOLAPt[j].coords[0];
      float lat = trackPts[i].LOLAPt[j].coords[1];
      float rad = trackPts[i].LOLAPt[j].coords[2];
      int id = trackPts[i].LOLAPt[j].s;

      Vector2 DEM_lonlat(lon, lat);
      Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(DEM_lonlat);

      int x = (int)DEM_pix[0];
      int y = (int)DEM_pix[1];

      PixelGray<float>  DEMVal = interpDEM(x, y);

      if ((id == ID) && (trackPts[i].valid)){
        //demPts.push_back((float)DEMVal);
        demPts[index] = (float)DEMVal; 
      }             
    }
    index++;
  }

  return demPts;
}

void SaveDEMPoints(vector< vector<LOLAShot> > &trackPts, string DEMFilename, string filename)

{    
  for (int k = 0; k < trackPts.size(); k++){
    vector<float> demPts = GetTrackPtsFromDEM(trackPts[k], DEMFilename, 3);
    string prefixTrackFilename =  prefix_from_filename(filename);  
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    SaveVectorToFile(demPts, string(trackFilename));     
  }
}


void SaveReflectancePoints(vector< vector<LOLAShot> >  &allTracks, float scaleFactor, string filename)
{

  FILE *fp;


  for (int k = 0; k < allTracks.size(); k++ ){


    string prefixTrackFilename =  prefix_from_filename(filename);  
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");


    for (int i = 0; i < allTracks[k].size(); i++){
      if (allTracks[k][i].valid == 1){
        fprintf(fp, "%f\n", scaleFactor*allTracks[k][i].reflectance);
      }
      else{
        fprintf(fp, "-1\n");
      }
    }

    fclose(fp);
    delete trackFilename;
  }

}

//saves the image points corresponding to a detectNum
void SaveImagePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename)
{

  FILE *fp;

  for (int k = 0; k < allTracks.size(); k++ ){

    string prefixTrackFilename = prefix_from_filename(filename); 
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");

    for (int i = 0; i < allTracks[k].size()-1; i++){


      int found = 0;
      for (int j = 0; j < allTracks[k][i].LOLAPt.size(); j++){
        if ((allTracks[k][i].LOLAPt[j].s == detectNum) && (allTracks[k][i].valid == 1)){
          found = 1;
          fprintf(fp, "%f\n", allTracks[k][i].imgPt[j].val);
        }
      }
      if (found == 0){
        fprintf(fp, "-1\n");
      }

    }

    fclose(fp);
    delete trackFilename;
  }

}

//saves the image points corresponding to a detectNum
void SaveAltitudePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename)
{

  FILE *fp;

  for (int k = 0; k < allTracks.size(); k++ ){

    string prefixTrackFilename = prefix_from_filename(filename); 
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");

    for (int i = 0; i < allTracks[k].size(); i++){

      int found = 0;
      for (int j = 0; j < allTracks[k][i].LOLAPt.size(); j++){
        if ((allTracks[k][i].LOLAPt[j].s == detectNum) && (allTracks[k][i].valid == 1)){
          found = 1;
          fprintf(fp, "%f\n", allTracks[k][i].LOLAPt[j].coords[2]);
        }
      }
      if (found == 0){
        fprintf(fp, "-1\n");
      }

    }

    fclose(fp);
    delete trackFilename;
  }

}


void SaveGCPoints(vector<vector<LOLAShot> > trackPts,  std::vector<std::string> imgFiles,  std::vector<int> overlapIndices, 
                  vector<Vector<float, 6> > optimalTransfArray, vector<float> optimalErrorArray, string gcpFilename)
{

  int index = 0;

  //for all features in the LOLA data
  for (int t=0; t<trackPts.size(); t++){
    for (int s=0; s<trackPts[t].size(); s++){
      if (trackPts[t][s].featurePtLOLA==1){

	float x = trackPts[t][s].LOLAPt[2].coords[0];
	float y = trackPts[t][s].LOLAPt[2].coords[1]; 
	float z = trackPts[t][s].LOLAPt[2].coords[2];
	float sigma_x = 1;
	float sigma_y = 1; 
	float sigma_z = 1;
 
	stringstream ss;
	ss<<index;
	string this_gcpFilename = gcpFilename+"_"+ss.str()+".gcp";
    
	FILE *fp = fopen(this_gcpFilename.c_str(), "w");
	
	fprintf(fp, "%f %f %f %f %f %f\n", x, y, z, sigma_x, sigma_y, sigma_z);
	for (int k = 0; k < overlapIndices.size(); k++){
	  float i = (optimalTransfArray[k][0]*trackPts[t][s].imgPt[2].x + optimalTransfArray[k][1]*trackPts[t][s].imgPt[2].y + optimalTransfArray[k][2]);
	  float j = (optimalTransfArray[k][3]*trackPts[t][s].imgPt[2].x + optimalTransfArray[k][4]*trackPts[t][s].imgPt[2].y + optimalTransfArray[k][5]);
          string filenameNoPath = imgFiles[overlapIndices[k]];
          int lastSlashPos = filenameNoPath.find_last_of("/");
          if (lastSlashPos != -1){
	    filenameNoPath.erase(0, lastSlashPos+1);
	  }
          fprintf(fp, "%s %f %f\n", filenameNoPath.c_str(), i, j);
	}
  
	fclose(fp);
       
        index++;
      }
    }
  }
  
}





