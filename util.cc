// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#ifndef UTIL_H
#define UTIL_H

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
#include <vw/Math.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::camera;
using namespace std;


void FindAndReplace( std::string& tInput, std::string tFind, std::string tReplace ) 
{ 

  size_t uPos = 0; 
  size_t uFindLen = tFind.length(); 
  size_t uReplaceLen = tReplace.length();
  
  if( uFindLen != 0 ){
    for( ;(uPos = tInput.find( tFind, uPos )) != std::string::npos; ){
      tInput.replace( uPos, uFindLen, tReplace );
      uPos += uReplaceLen;
    }
  }

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

  //printf("minLon = %f, minLat = %f, maxLon = %f, maxLat = %f\n", minLon, minLat, maxLon, maxLat);
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
       
       printf("lidar corners = %f %f %f %f\n", currCorners[0], currCorners[1], currCorners[2], currCorners[3]); 
       printf("image corners = %f %f %f %f\n", corners[0], corners[1], corners[2], corners[3]);       

       if(  ((corners(0)>currCorners(0)) && (corners(0)<currCorners(1))) //minlon corners in interval of currCorners 
	  ||((corners(1)>currCorners(0)) && (corners(1)<currCorners(1))) //maxlon corners in interval of currCorners
          ||((currCorners(0)>corners(0)) && (currCorners(0)<corners(1)))
          ||((currCorners(1)>corners(0)) && (currCorners(1)<corners(1)))) 
            
       {
         lonOverlap = 1;
       }
       if(  ((corners(2)>currCorners(2)) && (corners(2)<currCorners(3))) //minlat corners in interval of currCorners
	  ||((corners(3)>currCorners(2)) && (corners(3)<currCorners(3))) //maxlat corners in interval of currCorners
          ||((currCorners(2)>corners(2)) && (currCorners(2)<corners(3))) //minlat corners in interval of currCorners
	  ||((currCorners(3)>corners(2)) && (currCorners(3)<corners(3)))) //maxlat corners in interval of currCorners
	
       {
         latOverlap = 1;
       }
    
       cout<<"lonOverlap="<<lonOverlap<<", latOverlap="<<latOverlap<<endl; 
       cout<<"-----------------------------------------"<<endl;

       if ((lonOverlap == 1) && (latOverlap == 1)){
           overlapIndices.push_back(i);
       }
  }

  //cout<<overlapIndices<<endl;
  return overlapIndices;
}

#endif

