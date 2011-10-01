// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__
//
#ifndef UTIL_H
#define UTIL_H

#include <iostream>

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
int  ReadOverlapList(string lidarFilename, std::vector<int> &overlapIndices);

//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
Vector4 ComputeGeoBoundary(string cubFilename);


void ReadFileList(string camFileListFilename, vector<string> &camFileArray);

//determine if inputFiles is a text file containing a list of input files, one input file or a set of input files
//by reading the inputFiles extension. A text file containing a filename list *must* have extension .txt 
vector<string> AccessDataFilesFromInput(vector<string> &inputFiles);

//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> makeOverlapList(std::vector<std::string> inputFiles, Vector4 currCorners);
std::vector<int> makeOverlapListFromGeoTiff(std::vector<std::string> inputFiles, Vector4 currCorners);

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
SaveDEMErrors( const string&          filename, 
             const vector<Vector3>& locations, 
             const valarray<float>& errors,
             const vector<string>&  titles= vector<string>(),
             const string&          separator = ",",
             const string&          commentor = "#" );
void SaveStatistics (const string& filename, const vector<float>& errors, const vector<float> &histBins);
void SaveStatistics (const string& filename, const valarray<float>& errors, const vector<float> &histBins);

void ReadStatistics (const string& filename, vector<int>& hist, 
                     float *minError, float *maxError, float *avgError, int *numValidPts);
#endif

