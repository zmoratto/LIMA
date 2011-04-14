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

typedef struct CoregistrationParams{
  int matchingMode; //LIMA, LIDEM
  int reflectanceType;//NO, LAMBERT, LUNAR-LAMBERT
  int analyseFlag; 
  int useLOLAFeatures;
  int topPercentFeatures;
  int maxNumIter;
  int maxNumStarts;
  int displayResults;
  double noDataVal;
  float minConvThresh;
};


#endif /* COREGISTER_H */
