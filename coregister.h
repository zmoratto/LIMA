// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef COREGISTER_H
#define COREGISTER_H

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
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
#define LIMA 1
#define LIDEM 2

/*
struct ShotTime
{
  int year;
  int month;
  int day;
  int hour; //24hr format
  int min;
  float sec;
};

struct pointCloud
{
  int year;
  int month;
  int day;
  int hour; //24hr format
  int min;
  float sec;
  int s; //detector id
  Vector3 coords;
};

struct imgPoint
{
  double x;
  double y;
  double val;
};


// prototype for GetPointFromIndex
pointCloud GetPointFromIndex( vector<pointCloud> const & LOLAPts, int index);

struct LOLAShot
{
  int valid;
  vector<pointCloud> LOLAPt;
  vector<imgPoint> imgPt;
  float reflectance;
  float synthImage;
};

vector<float> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename, int ID);
vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID);
vector<Vector3> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename);
*/

#endif /* COREGISTER_H */
