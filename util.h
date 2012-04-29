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

void FindAndReplace( std::string& tInput, std::string tFind, std::string tReplace );
std::string GetFilenameNoExt(std::string const& filename);
std::string GetFilenameNoPath(std::string const& filename);
std::string GetFilenameExt(std::string const& filename);

void PrintOverlapList(std::vector<int>  overlapIndices);
void SaveOverlapList(string lidarFilename, std::vector<int> &overlapIndices);
void SaveOverlapList(string filename, std::vector<std::string> &filenames);
int  ReadOverlapList(string lidarFilename, std::vector<int> &overlapIndices);

//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
vw::Vector4 ComputeGeoBoundary( const std::string& );
vw::BBox2 ComputeGeoBBox( const std::string& );
vw::BBox2 ComputeGeoBBoxISIS( const std::string&, const vw::cartography::GeoReference& );
vw::BBox2 ComputeGeoBBoxISIS( const std::string& filename, const std::string& datumname );


std::vector<std::string> ReadFileList( const std::string& );
void ReadFileList( const std::string&, std::vector<std::string>& );

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

