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
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
#include "match.h"
#include "coregister.h"
#include "display.h"

vector<vector<LOLAShot> > CSVFileRead(string CSVFilename)
{
  string line;
  ifstream myfile (CSVFilename.c_str());
  int lineIndex = 0;
  //int shotIndex = 0;
 
  int trackIndex;
  vector<vector<LOLAShot> >trackPts;
  LOLAShot shot;
 
  pointCloud currPt;
  pointCloud prevPt;

  if (myfile.is_open())
  {
    while (! myfile.eof() )
    {
      getline (myfile,line);
      if(!myfile.eof()){ 
      if (lineIndex > 0){//skip header
     
        //float lon, lat, rad;
        char *temp = new char[160]; 
        char *lon = new char[160]; 
        char *lat = new char[160]; 
        char *rad = new char[160];
        char *detID = new char[5];
        char *tmpf = new char[200];
	
        //printf("line = %s\n", line.c_str()); 
	sscanf(line.c_str(), "%s %s %s %s %s %s %s %s %s %s  %s", 
	       temp, lon, lat, rad, tmpf, tmpf, tmpf, tmpf, tmpf, tmpf, detID);

        string time = temp;
        char year[5]; 
        char month[3]; 
        char day[3];
        char hour[3];
        char min[3];
        char sec[12];
        char s[2];
	
	string detIDs = detID;
     
        //NOTE: all of the following atoi where originally atof but where changed to remove compiler warnings - dtj, 2010_07_21     
        if (time.length() > 1){
	  size_t length;
	  length = time.copy(year, 4, 0);
          //printf("length = %d\n", length);
	  year[length] = '\0';
          currPt.year = atoi(year); 
       
	  length = time.copy(month, 2, 5);
          //printf("length = %d\n", length);
	  month[length] = '\0';
          currPt.month = atoi(month); 

	  length = time.copy(day, 2, 8);
	  day[length] = '\0';  
          currPt.day = atoi(day);

	  length = time.copy(hour, 2, 11);
	  hour[length] = '\0';
          currPt.hour = atoi(hour);

	  length = time.copy(min, 2, 14);
	  min[length] = '\0';
	  length = time.copy(sec, 11, 17);
	  sec[length] = '\0';
          currPt.sec = atof(sec);

          length = detIDs.copy(s, 1, 2);
          s[length] = '\0';
          currPt.s = atoi(s);
	  //printf("%s %s %s %s %s %s detID = %s\n", year, month, day, hour, min, sec, s);
	}
	
        Vector3 coords;         
        currPt.coords(0) = atof(lon);
        currPt.coords(1) = atof(lat);
	currPt.coords(2) = atof(rad);
        //cout<<"LOLAPt: "<<currPt.coords<<endl;
        
        if ((currPt.coords(0)!=0.0) && (currPt.coords(1)!=0.0) ){ //valid lidar point
     
	  if (lineIndex == 1){ //initialize first track
	      trackIndex = 0;
	      trackPts.resize(trackIndex+1);
             // printf("lineIndex = %d\n", lineIndex);
          }
          else{
            
	     if (GetTimeDiff(prevPt, currPt, 3000)){ //new track
	         trackPts[trackIndex].push_back(shot);//add last shot to the previous track
                 shot.LOLAPt.clear();
                 trackIndex++;
 	         trackPts.resize(trackIndex+1); //start new track
                 shot.LOLAPt.push_back(currPt);
	     }
	     else{ //same track
                 if (GetTimeDiff(prevPt, currPt, 0)){//new shot
	             trackPts[trackIndex].push_back(shot);
                     shot.LOLAPt.clear();
                     shot.LOLAPt.push_back(currPt);
                 }
	         else{ //same shot
                    shot.LOLAPt.push_back(currPt);
                 }
             }
	   }

      
           //copy current pc into prevPt
           prevPt.coords(0) = currPt.coords(0);
           prevPt.coords(1) = currPt.coords(1);
           prevPt.coords(2) = currPt.coords(2);
           prevPt.year = currPt.year;
           prevPt.month = currPt.month;
           prevPt.day = currPt.day;
           prevPt.hour = currPt.hour;
           prevPt.min = currPt.min;
           prevPt.sec = currPt.sec;   
           prevPt.s = currPt.s;
	}

        delete temp;
        delete lon;
	delete lat;
	delete rad;
	delete detID;
      } 
      lineIndex++; 
    }
    }
    myfile.close();
  }
  else
  {
  vw_throw( vw::IOErr() << "Unable to open track file \"" << CSVFilename << "\"" );
  }	
  
  return trackPts; 
}
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
  
  for(unsigned int k = 0; k < trackPts.size();k++){
    for(unsigned int i = 0; i < trackPts[k].size(); i++){
      
      trackPts[k][i].valid = 1; 

      LOLAPts = trackPts[k][i].LOLAPt;
      trackPts[k][i].imgPt.resize(LOLAPts.size());

      for (unsigned int j = 0; j < LOLAPts.size(); j++){
          
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
  for (unsigned int i = 0; i < trackPts.size(); i++){
    for (unsigned int k = 0; k < trackPts[i].LOLAPt.size(); k++){

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

  for(unsigned int k = 0; k < trackPts.size();k++){
    for(unsigned int i = 0; i < trackPts[k].size(); i++){
      if ((trackPts[k][i].valid == 1) && (trackPts[k][i].reflectance != 0) &&(trackPts[k][i].reflectance != -1)){//valid track and non-zero reflectance

        //update the nominator for the center point
        for (unsigned int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){
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


  for (unsigned int k = 0; k < allTracks.size(); k++ ){
    for (unsigned int i = 0; i < allTracks[k].size(); i++){
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
  for(unsigned int i = 0;i < LOLAPts.size(); i++){
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
  for (unsigned int i = 0; i < trackPts.size(); i++){
    demPts[index] = -1;
    for(unsigned int j = 0; j < trackPts[i].LOLAPt.size(); j++){

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
  for (unsigned int k = 0; k < trackPts.size(); k++){
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


  for (unsigned int k = 0; k < allTracks.size(); k++ ){


    string prefixTrackFilename =  prefix_from_filename(filename);  
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");


    for (unsigned int i = 0; i < allTracks[k].size(); i++){
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

  for (unsigned int k = 0; k < allTracks.size(); k++ ){

    string prefixTrackFilename = prefix_from_filename(filename); 
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");

    for (unsigned int i = 0; i < allTracks[k].size()-1; i++){


      int found = 0;
      for (unsigned int j = 0; j < allTracks[k][i].LOLAPt.size(); j++){
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

  for (unsigned int k = 0; k < allTracks.size(); k++ ){

    string prefixTrackFilename = prefix_from_filename(filename); 
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");

    for (unsigned int i = 0; i < allTracks[k].size(); i++){

      int found = 0;
      for (unsigned int j = 0; j < allTracks[k][i].LOLAPt.size(); j++){
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



void SaveGCPoints(vector<gcp> gcpArray,  string gcpFilename)
{

  for (unsigned int i = 0; i < gcpArray.size(); i++){
       stringstream ss;
       ss<<i;
       string this_gcpFilename = gcpFilename+"_"+ss.str()+".gcp";
    
       FILE *fp = fopen(this_gcpFilename.c_str(), "w");
	
       fprintf(fp, "%f %f %f %f %f %f\n", 
                    gcpArray[i].lon, gcpArray[i].lat, gcpArray[i].rad, 
                    gcpArray[i].sigma_lon, gcpArray[i].sigma_lat, gcpArray[i].sigma_rad);
       int numFiles = gcpArray[i].filename.size();
     
       for (int j = 0; j < numFiles-1; j++){
	 fprintf(fp,"%s %f %f\n", 
		    (char*)(gcpArray[i].filename[j].c_str()), gcpArray[i].x[j], gcpArray[i].y[j]);
       }
       if (numFiles > 0){
	 fprintf(fp, "%s %f %f", 
		 (char*)(gcpArray[i].filename[numFiles-1].c_str()), gcpArray[i].x[numFiles-1], gcpArray[i].y[numFiles-1]);
       }
       fclose(fp);
  }
   

}

//computes the min and max lon and lat of the LOLA data
//this function can be moved to a new util.c file
Vector4 FindMinMaxLat(vector<vector<LOLAShot> >trackPts)
{
  float minLat = 180;
  float maxLat = -180;
  float minLon = 180;
  float maxLon = -180;

  for (unsigned int i = 0; i < trackPts.size(); i++){
    for (unsigned int j = 0; j < trackPts[i].size(); j++){
      for(unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){
        float lon = trackPts[i][j].LOLAPt[k].coords[0];
        float lat = trackPts[i][j].LOLAPt[k].coords[1];

        if (lat < minLat){
          minLat = lat;
        }

        if (lon < minLon){
          minLon = lon;
        }

        if (lat > maxLat){
          maxLat = lat;
        }

        if (lon > maxLon){
          maxLon = lon;
        }
      }
    }
  }

  Vector4 coords;
  coords(0) = minLat; 
  coords(1) = maxLat; 
  coords(2) = minLon; 
  coords(3) = maxLon;

  return coords;

}

int GetTimeDiff(pointCloud prevPt, pointCloud currPt, float timeThresh)
{


  double prevTime = prevPt.hour*3600.0 + prevPt.min*60.0 + prevPt.sec;
  double currTime = currPt.hour*3600.0 + currPt.min*60.0 + currPt.sec;
 
  if (prevPt.year != currPt.year){//new orbit
      return (1);
  }  

  if (prevPt.month != currPt.month){ //new orbit
      return (1);
  } 
  
  if (prevPt.day != currPt.day){//new orbit
      return (1);
  } 
 
  if (currTime - prevTime > timeThresh){ //new orbit
     return (1);
  }
  else{
     return (0);
  }
}

//=============================================================================================================================
//this function can be removed - START
void SaveGCPoints(vector<vector<LOLAShot> > trackPts,  std::vector<std::string> imgFiles,  std::vector<int> overlapIndices, 
                  vector<Vector<float, 6> > optimalTransfArray, vector<float> optimalErrorArray, string gcpFilename)
{

  int index = 0;

  //for all features in the LOLA data
  for (unsigned int t=0; t<trackPts.size(); t++){
    for (unsigned int s=0; s<trackPts[t].size(); s++){
      if (trackPts[t][s].featurePtLOLA==1){

	float lon = trackPts[t][s].LOLAPt[2].coords[0];
	float lat = trackPts[t][s].LOLAPt[2].coords[1]; 
	float rad = trackPts[t][s].LOLAPt[2].coords[2]*1000;
	float sigma_x = 1;
	float sigma_y = 1; 
	float sigma_z = 1;
 
        //TO DO: compute the original imgPts from xyz - START
	Vector3 lon_lat_rad (lon,lat,rad);
	Vector3 xyz = lon_lat_radius_to_xyz(lon_lat_rad);
	
        //TO DO: compute the original imgPts from xyz - END

	stringstream ss;
	ss<<index;
	string this_gcpFilename = gcpFilename+"_"+ss.str()+".gcp";
    
	FILE *fp = fopen(this_gcpFilename.c_str(), "w");
	
	fprintf(fp, "%f %f %f %f %f %f\n", lon, lat, rad, sigma_x, sigma_y, sigma_z);
	
        for (unsigned int k = 0; k < overlapIndices.size(); k++){
          
          //boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
	  //double nodata_value = rsrc->nodata_read();
	  //DiskImageView<PixelGray<float> > isis_view( rsrc );
	  //int width = isis_view.cols();
	  //int height = isis_view.rows();

          string cubFilename =  imgFiles[overlapIndices[k]];
	  camera::IsisCameraModel model(cubFilename);
          Vector2 cub_pix = model.point_to_pixel(xyz);
	  float x = cub_pix[0];
	  float y = cub_pix[1];
          
          float i = (optimalTransfArray[k][0]*x + optimalTransfArray[k][1]*y + optimalTransfArray[k][2]);
	  float j = (optimalTransfArray[k][3]*x + optimalTransfArray[k][4]*y + optimalTransfArray[k][5]);

	  //float i = (optimalTransfArray[k][0]*trackPts[t][s].imgPt[2].x + optimalTransfArray[k][1]*trackPts[t][s].imgPt[2].y + optimalTransfArray[k][2]);
	  //float j = (optimalTransfArray[k][3]*trackPts[t][s].imgPt[2].x + optimalTransfArray[k][4]*trackPts[t][s].imgPt[2].y + optimalTransfArray[k][5]);
          
	  //print x and y vals before 
          cout<<"before"<<x<<" "<<y<<endl;
          //print x and y vals after
          cout<<"after"<<i<<" "<<j<<endl;

          string filenameNoPath = imgFiles[overlapIndices[k]];
          int lastSlashPos = filenameNoPath.find_last_of("/");
          if (lastSlashPos != -1){
	    filenameNoPath.erase(0, lastSlashPos+1);
	  }
          if (k == overlapIndices.size()-1){
	    fprintf(fp, "%s %f %f", filenameNoPath.c_str(), i, j);
          }
	  else{
	    fprintf(fp, "%s %f %f\n", filenameNoPath.c_str(), i, j);
	  }
	}
  
	fclose(fp);
       
        index++;
      }
    }
  }
  
}
//this function can be removed - END
