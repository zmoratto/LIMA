// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef PYR_TILING_H
#define PYR_TILING_H

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

typedef struct TilingParams
{
  int back_xl; //float coordinates to allow for subpixel precision 
  int back_yt;
  int back_xr;
  int back_yb;
  int horTileIndex;
  int verTileIndex;
  float backUpsampleFactor;
  float foreUpsampleFactor;
  string filename;
  string pcFilename;
  string accFilename;
};



void ComputePyramidTilingParams(int backImgWidth, int backImgHeight,  int tileSize,  
                                int &numTiles, Vector2 &offsetPix, int  &numPyrLevels);
void ComputeSubPyramidTilingParams(GeoReference const &foreGeo, GeoReference const &backGeo, Vector4 lonlatBB, int tileSize, Vector2 offsetPix,   Vector4 paddingParams,  int imageType,
				   std::vector<struct TilingParams> &tileParamsArray);
void ComputeTileParams(int orig_backImgWidth, int orig_backImgHeight, 
		       GeoReference const &foreGeo, GeoReference const &backGeo, 
		       int tileSize, Vector4 paddingParams, float foreNoDataVal,
		       int imageType, Vector4 lonlatBB, std::vector<struct TilingParams> &tileParamsArray);

#endif


