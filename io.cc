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
using namespace vw::camera;
using namespace std;

#include <math.h>
#include "match.h"
#include "coregister.h"
#include "tracks.h"
#include "io.h"



int ReadConfigFile(char *config_filename, struct CoregistrationParams *settings)
{
  int MAX_LENGTH = 5000;
  char line[MAX_LENGTH];
  ifstream configFile(config_filename);

  if (configFile.is_open()){
    printf("CONFIG FILE FOUND\n");
    
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MATCHING_MODE %d\n", &(settings->matchingMode));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "REFLECTANCE_TYPE %d\n", &(settings->reflectanceType));
   
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "ANALYSE_FLAG %d\n", &(settings->analyseFlag));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "USE_LOLA_FEATURES %d\n", &(settings->useLOLAFeatures));
   
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "TOP_PERCENT_FEATURES %d\n", &(settings->topPercentFeatures));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MAX_NUM_ITER %d\n", &(settings->maxNumIter));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MAX_NUM_STARTS %d\n", &(settings->maxNumStarts)); 

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "DISPLAY_RESULTS %d\n", &(settings->displayResults)); 

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "NO_DATA_VAL %lf\n", &(settings->noDataVal)); 

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "CONV_THRESH %f\n", &(settings->minConvThresh)); 
    
    configFile.close();

    return(1);
  }
  else{
    printf("configFile NOT FOUND\n");

    settings->matchingMode = LIMA;//LIDEM
    settings->reflectanceType = LUNAR_LAMBERT;// NO_REFL;//LAMBERT;
    settings->analyseFlag = 0;
    settings->useLOLAFeatures = 0;
    settings->topPercentFeatures = 10;
    settings->maxNumIter = 3;
    settings->maxNumStarts = 160;
    settings->displayResults = 0; 
    settings->noDataVal = -10000;
    settings->minConvThresh = 0.01;

    return(0);
  }
}

void PrintGlobalParams(struct CoregistrationParams *settings)
{

  printf("MATCHING_MODE %d\n", settings->matchingMode);
  printf("REFLECTANCE_TYPE %d\n", settings->reflectanceType);
  printf("ANALYSE_FLAG %d\n", settings->analyseFlag);
  printf("USE_LOLA_FEATURES %d\n", settings->useLOLAFeatures);
  printf("TOP_PERCENT_FEATURES %d\n", settings->topPercentFeatures);
  printf("MAX_NUM_ITER  %d\n", settings->maxNumIter);
  printf("MAX_NUM_STARTS %d\n", settings->maxNumStarts);
  printf("DISPLAY_RESULTS %d\n", settings->displayResults);
  printf("NO_DATA_VAL %f\n", settings->noDataVal);
  printf("MIN_CONV_THRESH %f\n", settings->minConvThresh);
}
//this function will be removed
int ReadModelParamsFile(string modelParamsFilename, struct ModelParams *params)
{

 int MAX_LENGTH = 5000;
 char line[MAX_LENGTH];
 ifstream modelParamsFile((char*)modelParamsFilename.c_str());

  if (modelParamsFile.is_open()){
    printf("MODEL PARAMS FILE FOUND\n");
    
    modelParamsFile.getline(line, MAX_LENGTH);

    modelParamsFile.getline(line, MAX_LENGTH);
    sscanf(line, "LIGHT_POS %lf %lf %lf\n", &(params->sunPosition[0]), &(params->sunPosition[1]), &(params->sunPosition[2])  );

    modelParamsFile.getline(line, MAX_LENGTH);
    sscanf(line, "VIEWER_POS %lf %lf %lf\n", &(params->spacecraftPosition[0]), &(params->spacecraftPosition[1]), &(params->spacecraftPosition[2]));

    params->exposureTime = 1.0;
    //params->rescalingParams[0] = 1;
    //params->rescalingParams[1] = 0;

    params->sunPosition[0] = 1000*(params->sunPosition[0]);
    params->sunPosition[1] = 1000*(params->sunPosition[1]);
    params->sunPosition[2] = 1000*(params->sunPosition[2]);

    params->spacecraftPosition[0] = 1000*(params->spacecraftPosition[0]);
    params->spacecraftPosition[1] = 1000*(params->spacecraftPosition[1]);
    params->spacecraftPosition[2] = 1000*(params->spacecraftPosition[2]);
    return (0);
  }

  return(1);
}

void PrintModelParams(struct ModelParams *params)
{
  cout<<"LIGHT_POS[0]"<<params->sunPosition[0]<<endl;
  cout<<"LIGHT_POS[1]"<<params->sunPosition[1]<<endl;
  cout<<"LIGHT_POS[2]"<<params->sunPosition[2]<<endl;

  cout<<"VIEWER_POS[0]"<<params->spacecraftPosition[0]<<endl;
  cout<<"VIEWER_POS[1]"<<params->spacecraftPosition[1]<<endl;
  cout<<"VIEWER_POS[2]"<<params->spacecraftPosition[2]<<endl;
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

void SaveVectorToFile(vector<float> v, string filename)
{
 

  FILE *fp;

  fp = fopen(filename.c_str(), "w");

  for (int i = 0; i < v.size()-1; i++){
    fprintf(fp, "%f\n", v[i]);
  }
  fprintf(fp, "%f", v[v.size()-1]);
  fclose(fp);
}


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
              printf("lineIndex = %d\n", lineIndex);
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

  else cout << "Unable to open file";

  
  return trackPts; 
}

//computes the min and max lon and lat of the LOLA data
//this function can be moved to a new util.c file
Vector4 FindMinMaxLat(vector<vector<LOLAShot> >trackPts)
{
  float minLat = 180;
  float maxLat = -180;
  float minLon = 180;
  float maxLon = -180;

  for (int i = 0; i < trackPts.size(); i++){
    for (int j = 0; j < trackPts[i].size(); j++){
      for(int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){
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

void printOverlapList(std::vector<int>  overlapIndices)
{
  printf("numOverlapping images = %d\n", (int)(overlapIndices.size()));
    for (int i = 0; i < overlapIndices.size(); i++){
      printf("%d ", overlapIndices[i]);
    }
}

//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
Vector4 ComputeGeoBoundary(string cubFilename)
{
  GeoReference moonref( Datum("D_MOON"), identity_matrix<3>() );
  boost::shared_ptr<IsisCameraModel> isiscam( new IsisCameraModel( cubFilename ) );
  BBox2 camera_boundary =
    camera_bbox( moonref, boost::shared_dynamic_cast<CameraModel>(isiscam),
		 isiscam->samples(), isiscam->lines() );
 
  Vector4 corners;
 

  float minLon = camera_boundary.min()[0];
  float minLat = camera_boundary.min()[1];
  float maxLon = camera_boundary.max()[0];
  float maxLat = camera_boundary.max()[1];

  printf("minLon = %f, minLat = %f, maxLon = %f, maxLat = %f\n", minLon, minLat, maxLon, maxLat);
  if (maxLat<minLat){
      float temp = minLat;
      minLat = maxLat;
      maxLat = temp;    
  }

  if (maxLon<minLon){
     float temp = minLon;
     minLon = maxLon;
     maxLon = temp;    
  }

  corners(0) = minLon;
  corners(1) = maxLon;
  corners(2) = minLat;
  corners(3) = maxLat;

  return corners;
}
//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> makeOverlapList(std::vector<std::string> inputFiles, Vector4 currCorners) {
  
  std::vector<int> overlapIndices;
 
  for (unsigned int i = 0; i < inputFiles.size(); i++){

       int lonOverlap = 0;
       int latOverlap = 0; 
     
       Vector4 corners = ComputeGeoBoundary(inputFiles[i]);

       printf("corners = %f %f %f %f\n", corners[0], corners[1], corners[2], corners[3]);       

       if(  ((corners(0)>currCorners(0)) && (corners(0)<currCorners(1))) //minlon in interval 
	    ||((corners(1)>currCorners(0)) && (corners(1)<currCorners(1)))) //maxlon in interval
       {
         lonOverlap = 1;
       }
       if(  ((corners(2)>currCorners(2)) && (corners(2)<currCorners(3))) //minlat in interval 
	    ||((corners(3)>currCorners(2)) && (corners(3)<currCorners(3)))) //maxlat in interval
       {
         latOverlap = 1;
       }
     
       if ((lonOverlap == 1) && (latOverlap == 1)){
           overlapIndices.push_back(i);
       }
  }

  return overlapIndices;
}












