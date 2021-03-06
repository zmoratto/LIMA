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
#include <boost/filesystem.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include "util.h"

namespace fs = boost::filesystem;

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::camera;
using namespace std;


void PrintOverlapList( const std::vector<int>& overlapIndices ) {
  cout << "numOverlapping = " << overlapIndices.size() << endl;
  std::vector<int>::const_iterator i;
  for( i = overlapIndices.begin(); i != overlapIndices.end(); ++i){
    cout << "index: " << *i << endl;
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

BBox2 ComputeGeoBBoxISIS( const string& f, const GeoReference& g ) {
  boost::shared_ptr<IsisCameraModel> isiscam( new IsisCameraModel( f ) );
  BBox2 bbox = camera_bbox( g, boost::shared_dynamic_cast<CameraModel>(isiscam),
		                    isiscam->samples(), isiscam->lines() ); 
  return bbox;
}

BBox2 ComputeGeoBBoxISIS( const string& filename, const string& datumname ) {
  GeoReference georef( Datum( datumname ), identity_matrix<3>() );
  return ComputeGeoBBoxISIS( filename, georef );
}

BBox2 ComputeGeoBBox( const string& f ) {
  GeoReference g;
  vw_log().console_log().rule_set().add_rule(-1,"fileio"); // Silence GDAL
  bool r = read_georeference( g, f );
  vw_settings().reload_config();                           // unsilence it
  if( !r ){
    // Can't read a GeoReference from the file, so we will now assume that
    // it is an unprojected ISIS cube.  We also assume the Moon, which is perhaps
    // not the best.
    return ComputeGeoBBoxISIS( f, "D_MOON" );
  }
  
  DiskImageView<PixelGray<float> > image( f );
  BBox2 image_bbox = g.pixel_to_lonlat_bbox( bounding_box(image) );
  return image_bbox;
}

Vector4 ComputeGeoBoundary( const string& cubFilename ) {

  BBox2 camera_boundary = ComputeGeoBBox( cubFilename );
  Vector4 corners;
  corners(0) = camera_boundary.min().x();
  corners(1) = camera_boundary.max().x();
  corners(2) = camera_boundary.min().y();
  corners(3) = camera_boundary.max().y();
  return corners;
}

//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> makeOverlapList( const std::vector<std::string>& inputFiles,
                                  const vw::BBox2&                bbox ) {
  std::vector<int> overlapIndices;
  for( unsigned int i = 0; i < inputFiles.size(); ++i ){
    BBox2 image_bbox = ComputeGeoBBox( inputFiles[i] );
    if( bbox.intersects(image_bbox) ){ overlapIndices.push_back(i); }
  }
  return overlapIndices;
} 
    

std::vector<int> makeOverlapList( const std::vector<std::string>& inputFiles, 
                                  const Vector4&                  currCorners) {
  Vector2 min( currCorners(0), currCorners(2) );
  Vector2 max( currCorners(1), currCorners(3) );
  BBox2 bbox( min, max );
  return makeOverlapList( inputFiles, bbox );
}

// Since makeOverlapList() is now completely generalized, this isn't needed anymore,
// but is retained for historical interface reasons.
std::vector<int> makeOverlapListFromGeoTiff( const vector<string>& inputFiles, 
                                             const Vector4&        currCorners) {
  return makeOverlapList( inputFiles, currCorners );
}

void SaveOverlapList( const string& filename, const std::vector<int>& overlapIndices ) {
  WriteVectorTo( filename, overlapIndices );
}

void SaveOverlapList( const string& filename, const std::vector<std::string>& filenames ) {
  WriteVectorTo( filename, filenames );
}

int ReadOverlapList( const string& filename, std::vector<int>& overlapIndices ) {
  try{
    overlapIndices = ReadVectorFrom<int>( filename );
    return 1;
  }
  catch( const vw::ArgumentErr& error ){
    return 0;
  }
}


void SaveDEMErrors( const string&          filename, 
                    const vector<Vector3>& locations, 
                    const valarray<float>& errors,
                    const vector<string>&  titles,
                    const string&          separator,
                    const string&          commentor) {
  ofstream file( filename.c_str() );
   
  if( !file ) {
    vw_throw( ArgumentErr() << "Can't open error output file \"" << filename << "\"" );
  }
   
  if( locations.size() != errors.size() ) {
    vw_throw( ArgumentErr() << "The there are a different number of locations (" 
                            << locations.size() << ") than errors (" 
                            << errors.size() <<") which is a problem." );
  }
   
  if( titles.size() > 0 ){
    file << commentor << " ";
    string title = boost::join( titles, separator );
    file << title << endl;
  }
   
  for( unsigned int i = 0; i < locations.size(); i++ ){
    file << locations[i].x() << separator
         << locations[i].y() << separator
         << fixed
         << locations[i].z() << separator
         << errors[i] << endl;
  }
  file.close();
}


void SaveStatistics( const string&        filename, 
                     const vector<float>& errors, 
                     const vector<float>& histBins) {
  valarray<float> valerrors(errors.size());
  for( unsigned int i = 0; i < errors.size(); ++i ){
    valerrors[i] = errors[i];
  }
  SaveStatistics( filename, valerrors, histBins );
}

void SaveStatistics( const string&          filename, 
                     const valarray<float>& errors, 
                     const vector<float>&   histBins) {
  float minError = 10000000.0; 
  float maxError = -10000000.0;
  float avgError = 0.0;
  int numValidPts = 0;
  vector<float> errorHist( histBins.size() );

  for( unsigned int i = 0; i < errors.size(); ++i ){
    if( errors[i] < 0 ){
      cout << "writeStats: inavlid error" << errors[i] << endl;
      continue;
      }
	float stdv = errors[i];
	
    //update the number of valid points
    ++numValidPts;
        
    //update the bins
    for( unsigned int k = 1; k < histBins.size(); ++k ){
	  if( (stdv > histBins[k-1]) && (stdv <= histBins[k]) ){
	    ++errorHist[k-1];
        break;
      }
    }
	if( stdv > histBins.back() ){ errorHist.back()++; }

    //update the average
    avgError += stdv;

    //update the max error
	if( stdv > maxError ){ maxError = stdv; }  
        
    //update the min error
	if( stdv< minError ){ minError = stdv; }
  }  

  if( numValidPts ){ avgError /= numValidPts; }

  ofstream file( filename.c_str() );

  if( !file ) {
      vw_throw( ArgumentErr() << "Can't open statistics output file \"" << filename << "\"" );
  }

  file << "minError= "    << minError    << endl;
  file << "maxError= "    << maxError    << endl;
  file << "avgError= "    << avgError    << endl;
  file << "numValidPts= " << numValidPts << endl;

  for( unsigned int i = 0; i < histBins.size(); ++i){
    file << "hist_" << i << "= " << errorHist[i] << endl;
  }
  file.close();
}


void ReadStatistics( const string&      filename,
                           vector<int>& hist, 
                           float*       minError, 
                           float*       maxError, 
                           float*       avgError, 
                           int*         numValidPts) {
  ifstream file( filename.c_str() );

  if( !file ) {
      vw_throw( ArgumentErr() << "Can't open statistics input file \"" << filename << "\"" );
  }

  int val;
  string temp;
  float l_avgError;
  float l_minError;
  float l_maxError;
  int l_numValidPts;

  file >> temp >> l_minError;
  file >> temp >> l_maxError;
  file >> temp >> l_avgError;
  file >> temp >> l_numValidPts;
   
  while( file >> temp >> val ){ hist.push_back(val); }
  file.close();
  
  *minError = l_minError;
  *maxError = l_maxError;
  *avgError = l_avgError;
  *numValidPts = l_numValidPts;
}



void ReadFileList( const string& filepath, vector<string>& fileArray ) {
  fileArray.clear();
  fileArray = ReadVectorFrom<string>( filepath );
  return;
}

vector<string> AccessDataFilesFromInput(const string& inputFile) {
  const vector<string> v(1,inputFile);
  return AccessDataFilesFromInput( v );
}

vector<string> AccessDataFilesFromInput(const vector<string>& inputFiles) {
  vector<string> v;
  if( inputFiles.empty() ){ 
    vw_throw( ArgumentErr() << "There are no input files to read." );
  }
  for( vector<string>::const_iterator i = inputFiles.begin();
       i != inputFiles.end();
       ++i ){
    fs::path p( *i );
    if( p.extension() == string(".txt") ){
      vector<string> list = ReadVectorFrom<string>( p );
      v.insert( v.end(), list.begin(), list.end() );
    }
    else{ v.push_back( *i ); }
  }
  return v;
}
