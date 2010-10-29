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
#include "match.h"
#include "coregister.h"
#include "tracks.h"


typedef struct GlobalParams{
  int matchingMode; //LIMA, LIDEM
  int reflectanceType;//NO, LAMBERT, LUNAR-LAMBERT
  int analyseFlag; 
  int maxNumIter;
  int maxNumStarts;
};

int ReadConfigFile(char *config_filename, struct GlobalParams *settings)
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
    sscanf(line, "MAX_NUM_ITER %d\n", &(settings->maxNumIter));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MAX_NUM_STARTS %d\n", &(settings->maxNumStarts)); 

    //configFile.getline(line, MAX_LENGTH);
    
    configFile.close();

    return(1);
  }
  else{
    printf("configFile NOT FOUND\n");

    settings->matchingMode = LIMA;//LIDEM
    settings->reflectanceType = LUNAR_LAMBERT;// NO_REFL;//LAMBERT;
    settings->analyseFlag = 0;
    settings->maxNumIter = 3;
    settings->maxNumStarts = 160; 

    return(0);
  }
}

void PrintGlobalParams(struct GlobalParams *settings)
{

  printf("MATCHING_MODE %d\n", settings->matchingMode);
  printf("REFLECTANCE_TYPE %d\n", settings->reflectanceType);
  printf("ANALYSE_FLAG %d\n", settings->analyseFlag);
  printf("MAX_NUM_ITER  %d\n", settings->maxNumIter);
  printf("MAX_NUM_STARTS %d\n", settings->maxNumStarts);
 
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
    //NOTE: all of the following atoi where orignally atof but where changed to remove compiler warnings - dtj, 2010_07_21
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










