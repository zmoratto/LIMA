// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
//
#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <vw/Core/Exception.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>

#include <vw/Cartography.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

void ReadFileList( const std::string&, std::vector<std::string>& );
void PrintOverlapList( const std::vector<int>& );
void SaveOverlapList( const std::string&, const std::vector<int>& );
void SaveOverlapList( const std::string&, const std::vector<std::string>& );
int ReadOverlapList( const std::string&, std::vector<int>&);

//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
vw::Vector4 ComputeGeoBoundary( const std::string& );
vw::BBox2 ComputeGeoBBox( const std::string& );
vw::BBox2 ComputeGeoBBoxISIS( const std::string&, const vw::cartography::GeoReference& );
vw::BBox2 ComputeGeoBBoxISIS( const std::string& filename, const std::string& datumname );

//determine if inputFiles is a text file containing a list of input files, one input file or a set of input files
//by reading the inputFiles extension. A text file containing a filename list *must* have extension .txt 
std::vector<std::string> AccessDataFilesFromInput(const std::string& );
std::vector<std::string> AccessDataFilesFromInput(const std::vector<std::string>& );

//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> makeOverlapList( const std::vector<std::string>&, const vw::Vector4& );
std::vector<int> makeOverlapList( const std::vector<std::string>&, const vw::BBox2& );
std::vector<int> makeOverlapListFromGeoTiff( const std::vector<std::string>&, 
                                             const vw::Vector4& );

/* Stream manipulator to ignore until the end of line
 * Use like: std::cin >> ignoreLine;
 */
template <class charT, class traits>
inline
std::basic_istream<charT, traits>&
ignoreLine (std::basic_istream<charT,traits>& stream)
	{
	// skip until end of line
	stream.ignore( std::numeric_limits<int>::max(), stream.widen('\n') );

	return stream;
	}

/* Stream manipulator to ignore one character
 * Use like: std::cin >> ignoreOne;
 */
template <class charT, class traits>
inline
std::basic_istream<charT, traits>&
ignoreOne(std::basic_istream<charT,traits>& stream)
	{
	stream.ignore();

	return stream;
	}

inline Vector3 find_centroid( const vector<Vector3>& points )
  {
  Vector3 centroid;
  for (unsigned int i = 0; i < points.size(); i++)
    {
    centroid += points[i];
    }
  centroid /= points.size();
  return centroid;
  }

template <class T>
inline std::vector<T> ReadVectorFrom( std::istream& stream ){
  std::vector<T> v;
  T item;
  while( stream >> item ){ v.push_back( item ); }
  return v;
}

template <class T>
inline std::vector<T> ReadVectorFrom( const std::string& s ){
  std::ifstream file( s.c_str() );
    if (!file){
      vw_throw( vw::ArgumentErr() << "Can't open file " << s );
    }
  return ReadVectorFrom<T>( file );
}

template <class T>
inline std::vector<T> ReadVectorFrom( const boost::filesystem::path& p ){
  return ReadVectorFrom<T>( p.string() );
}

template <class T>
inline std::ostream& WriteVectorTo( std::ostream& stream, const std::vector<T>& v){
  if( v.empty() ){ return stream; }
  for( typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i ){
    stream << *i << std::endl;
  }
  return stream;
}

template <class T>
inline void WriteVectorTo( const std::string& s, const std::vector<T>& v){
  if( v.empty() ){
    vw_throw( vw::ArgumentErr() << "WriteVectorTo: the input vector is empty." );
  }

  std::ofstream file( s.c_str() );
  if( !file ) {
    vw_throw( vw::ArgumentErr() << "Can't open file \"" << s << "\"" );
  }

  WriteVectorTo( file, v );
  file.close();
  return;
}

template <class T>
inline void WriteVectorTo( const boost::filesystem::path& p, const std::vector<T>& v){
  WriteVectorTo( p.string(), v );
  return;
}



// Writes out the locations and errors to a file
void
SaveDEMErrors( const std::string&              filename, 
               const std::vector<vw::Vector3>& locations, 
               const std::valarray<float>&     errors,
               const std::vector<std::string>& titles= vector<string>(),
               const std::string&              separator = ",",
               const std::string&              commentor = "#" );

void SaveStatistics( const std::string&        filename, 
                     const std::vector<float>& errors, 
                     const std::vector<float>& histBins);

void SaveStatistics( const std::string&          filename, 
                     const std::valarray<float>& errors, 
                     const std::vector<float>&   histBins);

void ReadStatistics( const std::string&      filename, 
                           std::vector<int>& hist, 
                           float*            minError, 
                           float*            maxError, 
                           float*            avgError, 
                           int*              numValidPts );
#endif

