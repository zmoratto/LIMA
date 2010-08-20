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

//computes the scale factor for all tracks at once
float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts)
{
  float nominator = 0.0;
  int numValidPts = 0;
  float scaleFactor = 1;

   for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      if ((trackPts[k][i].valid == 1) && (trackPts[k][i].reflectance != 0)){//valid track and non-zero reflectance

        //update the nominator for the center point
	for (int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){
	  if (trackPts[k][i].LOLAPt[j].s == 3){ 
            printf("img = %f, refl = %f\n", trackPts[k][i].imgPt[j].val, trackPts[k][i].reflectance);
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

void ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks, ModelParams modelParams, GlobalParams globalParams)
{

  vector<pointCloud> LOLAPts;

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


vector<float> ComputeSyntImgPts(float scaleFactor, vector<vector<LOLAShot > >&trackPts)
{
  
 vector<float>synthImg;
 float thisSynthImg;
 int index = 0;
 int num_allImgPts = 0;
 for(int k = 0; k<trackPts.size();k++){
    num_allImgPts += trackPts[k].size();
  }
  synthImg.resize(num_allImgPts);  

 for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){

      if (trackPts[k][i].valid == 1){//valid track        
          synthImg[index] = scaleFactor*trackPts[k][i].reflectance; 
      }
      else{
	synthImg[index] = -1;
      }

      index++;
    } //i
 } //k

  return synthImg;
}

void SaveReflectance(vector< vector<LOLAShot> >  &allTracks, string filename)
{

  FILE *fp;
  fp = fopen(filename.c_str(), "w");
    for (int k = 0; k < allTracks.size(); k++ ){
      for (int i = 0; i < allTracks[k].size(); i++){
        if (allTracks[k][i].valid == 1){
	  fprintf(fp, "%f ", allTracks[k][i].reflectance);
	}
        else{
           fprintf(fp, "-1 ");
        }
      }
      fprintf(fp, "\n");
    }
  fclose(fp);
}

//saves the image points corresponding to a detectNum
void SaveImagePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename)
{

  FILE *fp;
  fp = fopen(filename.c_str(), "w");
    for (int k = 0; k < allTracks.size(); k++ ){
      for (int i = 0; i < allTracks[k].size(); i++){
        
        int found = 0;
	for (int j = 0; j < allTracks[i][k].LOLAPt.size(); j++){
	  if (allTracks[i][k].LOLAPt[j].s == detectNum){
	     found = 1;
	     fprintf(fp, "%f ", allTracks[k][i].imgPt[j].val);
	  }
	}
        if (found == 0){
           fprintf(fp, "-1");
        }

      }
      fprintf(fp, "\n");
    }
  fclose(fp);
}

