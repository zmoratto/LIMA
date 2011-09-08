// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/algorithm/string/join.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include "util.h"

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

///Returns the file extension
std::string GetFilenameExt(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(0, index+1);
  return result;
}

/// Erases a file suffix if one exists and returns the base string
std::string GetFilenameNoExt(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
      result.erase(index, result.size());
  return result;
}
/*
/// Erases a file suffix if one exists and returns the base string less3 characters
static std::string prefix_less3_from_filename(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
    result.erase(index-3, result.size()+3);
  return result;
}
*/

/// Erases a file path if one exists and returns the base string 
std::string GetFilenameNoPath(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind("/");
  if (index != -1)
    result.erase(0, index+1);
  return result;
}

void PrintOverlapList(std::vector<int>  overlapIndices)
{
  printf("numOverlapping images = %d\n", (int)(overlapIndices.size()));
    for (unsigned int i = 0; i < overlapIndices.size(); i++){
      printf("%d ", overlapIndices[i]);
    }
}

//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
Vector4 ComputeGeoTiffBoundary(string filename)
{

  Vector4 corners;
  
  GeoReference Geo;
  read_georeference(Geo, filename);
  DiskImageView<PixelGray<float> > DEM(filename);

  Vector2 pixel_tl(0,0);
  Vector2 lonlat_tl = Geo.pixel_to_lonlat(pixel_tl);	
  Vector2 pixel_br(DEM.cols()-1, DEM.rows()-1);
  Vector2 lonlat_br = Geo.pixel_to_lonlat(pixel_br);		
 
  float minLon = lonlat_tl[0];
  float minLat = lonlat_tl[1];
  float maxLon = lonlat_br[0];
  float maxLat = lonlat_br[1];

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
   
       //printf("lidar corners = %f %f %f %f\n", currCorners[0], currCorners[1], currCorners[2], currCorners[3]); 
       //printf("image corners = %f %f %f %f\n", corners[0], corners[1], corners[2], corners[3]);       

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
    
       //cout<<"lonOverlap="<<lonOverlap<<", latOverlap="<<latOverlap<<endl; 
       //cout<<"-----------------------------------------"<<endl;

       if ((lonOverlap == 1) && (latOverlap == 1)){
           overlapIndices.push_back(i);
       }
  }

  //cout<<overlapIndices<<endl;
  return overlapIndices;
}

//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> makeOverlapListFromGeoTiff(std::vector<std::string> inputFiles, Vector4 currCorners) 
{  
  std::vector<int> overlapIndices;
 
  for (unsigned int i = 0; i < inputFiles.size(); i++){

       int lonOverlap = 0;
       int latOverlap = 0; 
     
      
       Vector4 corners = ComputeGeoTiffBoundary(inputFiles[i]);

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

void SaveOverlapList(string filename, std::vector<int> &overlapIndices)
{
   ofstream file( filename.c_str() );
   //cout<<"numImgsOverlap="<<overlapIndices.size()<<endl;
   if (overlapIndices.size() > 0){
     for (int i = 0; i < overlapIndices.size()-1; i++){
       file<<overlapIndices[i]<<endl;
     }
   
     file<<overlapIndices[overlapIndices.size()-1];
   }
   else{
      file<<-1;
   }
   file.close();
   
}
int ReadOverlapList(string filename, std::vector<int> &overlapIndices)
{
   int fileFound = 0;
   overlapIndices.clear();
   ifstream file( filename.c_str() );
   if (!file){
     cout<<"file not found"<<endl;
     return fileFound;
   }
   else{
     cout<<"file found"<<endl;
     fileFound = 1;
     while (!file.eof()){
       int index;
       file>>index;  
       overlapIndices.push_back(index);
     }
     file.close();
     return fileFound;
   }
}

void SaveDEMErrors( const string& filename, 
		  const vector<Vector3>& locations, 
		  const valarray<float>& errors,
		  const vector<string>&  titles,
		  const string& separator,
		  const string& commentor)
{

   ofstream file( filename.c_str() );
   
   if( !file ) {
     vw_throw( ArgumentErr() << "Can't open error output file \"" << filename << "\"" );
   }
   
   if( locations.size() != errors.size() ) {
     vw_throw( ArgumentErr() 
	       << "The there are a different number of locations (" 
	       << locations.size() << ") than errors (" 
	       << errors.size() <<") which is a problem." );
   }
   
   if( titles.size() > 0 ){
     file << commentor << " ";
     string title = boost::join( titles, separator );
     file << title << endl;
   }
   
   for( unsigned int i = 0; i < locations.size(); i++ ){
     file  << locations[i].x() << separator
	   << locations[i].y() << separator
	   << fixed
	   << locations[i].z() << separator
	   << errors[i] << endl;
  }
   
   
   file.close();
}


void SaveStatistics (const string& filename, const vector<float>& errors, const vector<float> &histBins)
{
    float minError = 10000000.0; 
    float maxError = -10000000.0;
    float avgError = 0.0;
    int numValidPts = 0;

    int numBins = histBins.size();
    vector<float> errorHist;
    errorHist.resize(numBins);
  

    for ( unsigned int i = 0; i < errors.size(); i++ ){

      if (errors[i]<0){
	cout<<"writeStats: inavlid error"<<errors[i]<<endl;
      }
      float stdv;

      if (errors[i]>=0){
     
	float stdv = errors[i];
	
        //update the number of valid points
        numValidPts++;
        
        //update the bins
        for (int k = 0; k < histBins.size(); k++){
	  if ((stdv>histBins[k]) && (stdv<= histBins[k+1])){
	    errorHist[k]++;
           }
        }
	if (stdv > histBins.back()){
	  errorHist[histBins.size()-1]++;
	}

        //update the average
        avgError = avgError + stdv;

        //update the max error
	if (stdv> maxError){
	  maxError = stdv; 
	}  
        
        //update the min error
	if (stdv< minError){
	  minError = stdv; 
	}

      }  

    }

    if (numValidPts){
      avgError = avgError/numValidPts;
    }

    ofstream file( filename.c_str() );
    file<<"minError= "<<minError<<endl;
    file<<"maxError= "<<maxError<<endl;
    file<<"avgError= "<<avgError<<endl;
    file<<"numValidPts= "<<numValidPts<<endl;

    for (int i = 0; i<numBins; i++){
      file<<"hist_"<<i<<"= "<<errorHist[i]<<endl;
    }
  
    file.close();
}

void SaveStatistics (const string& filename, const valarray<float>& errors, const vector<float> &histBins)
{
    float minError = 10000000.0; 
    float maxError = -10000000.0;
    float avgError = 0.0;
    int numValidPts = 0;
 
    int numBins = histBins.size();
    vector<float> errorHist;
    errorHist.resize(numBins);

    for ( unsigned int i = 0; i < errors.size(); i++ ){

      if (errors[i]<0){
	cout<<"writeStats: inavlid error"<<errors[i]<<endl;
      }
      float stdv;

      if (errors[i]>=0){
     
	float stdv = errors[i];
	
        //update the number of valid points
        numValidPts++;
        
        //update the bins
        for (int k = 0; k < histBins.size(); k++){
	  if ((stdv>histBins[k]) && (stdv<= histBins[k+1])){
	    errorHist[k]++;
           }
        }
	if (stdv > histBins.back()){
	  errorHist[histBins.size()-1]++;
	}

        //update the average
        avgError = avgError + stdv;

        //update the max error
	if (stdv> maxError){
	  maxError = stdv; 
	}  
        
        //update the min error
	if (stdv< minError){
	  minError = stdv; 
	}

      }  

    }

    if (numValidPts){
      avgError = avgError/numValidPts;
    }

    ofstream file( filename.c_str() );
    file<<"minError= "<<minError<<endl;
    file<<"maxError= "<<maxError<<endl;
    file<<"avgError= "<<avgError<<endl;
    file<<"numValidPts= "<<numValidPts<<endl;

    for (int i = 0; i<numBins; i++){
      file<<"hist_"<<i<<"= "<<errorHist[i]<<endl;
    }
  
    file.close();
}


void ReadStatistics (const string& filename, vector<int>& hist, 
                     float *minError, float *maxError, float *avgError, int *numValidPts)
{
  ifstream file;
  file.open(filename.c_str());

  int val;
  string temp;
  float l_avgError;
  float l_minError;
  float l_maxError;
  int l_numValidPts;

  if (file.is_open()) {
    file >> temp;
    file >> l_minError;
    file >> temp;
    file >> l_maxError;
    file >> temp;
    file >> l_avgError;
    file >> temp;
    file >> l_numValidPts;
   
    while (!file.eof()) {
      file >> temp;
      file >> val;
      hist.push_back(val);
      //cout<<val<<endl;
    }
  
    *minError = l_minError;
    *maxError = l_maxError;
    *avgError = l_avgError;
    *numValidPts = l_numValidPts;
 
  }
  file.close();
    
}
// this should be moved to util
void ReadFileList(string fileListFilename, vector<string> &fileArray)
{

 ifstream file;
 file.open(fileListFilename.c_str());
 string filename;
 if (file.is_open()) {
    while (!file.eof()) {
      file >> filename;
      fileArray.push_back(filename);
      cout<<filename<<endl;
    }
  }
  file.close();

}
