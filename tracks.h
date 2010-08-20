// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef TRACKS_H
#define TRACKS_H

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

//vector<float> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename, int ID);
//vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID);
//vector<Vector3> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename);
//computes the scale factor for all tracks at once
float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts);
float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts);
void ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks, ModelParams modelParams, GlobalParams globalParams);
vector<float> ComputeSyntImgPts(float scaleFactor, vector<vector<LOLAShot > >&trackPts);
void SaveReflectance(vector< vector<LOLAShot> >  &allTracks, string filename);
void SaveImagePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename);

template <class ViewT>
void 
GetAllPtsFromImage(vector<vector<LOLAShot > > &trackPts,  ImageViewBase<ViewT> const& DRG, GeoReference const &DRGGeo)
{

  vector<pointCloud> LOLAPts;

 
  ImageViewRef<float> interpDRG;
  if ( IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(DRG.impl(),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );

    //cout << "NOT masked" <<endl;
  } else {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(apply_mask(DRG.impl()),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }

  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
   
      LOLAPts = trackPts[k][i].LOLAPt;
      trackPts[k][i].imgPt.resize(LOLAPts.size());

      for (int j = 0; j < LOLAPts.size(); j++){
          
	    float lon = LOLAPts[j].coords[0];
	    float lat = LOLAPts[j].coords[1];
	    float rad = LOLAPts[j].coords[2];
     
	    Vector2 DEM_lonlat(lon, lat);
	    Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);
	  
	    float x = DRG_pix[0];
	    float y = DRG_pix[1];

	    trackPts[k][i].imgPt[j].val = interpDRG(x, y);
	    trackPts[k][i].imgPt[j].x = DRG_pix[0];
	    trackPts[k][i].imgPt[j].y = DRG_pix[1];
	
      }

      trackPts[k][i].valid = 0; 

      //check for valid shot
      pointCloud centerPt  = GetPointFromIndex( LOLAPts, 3);
      pointCloud topPt     = GetPointFromIndex( LOLAPts, 2);
      pointCloud leftPt    = GetPointFromIndex( LOLAPts, 1);

      if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1) && (LOLAPts.size() <=5)){//valid LOLA shot
          
          trackPts[k][i].valid = 1;
	
      }
     
    }//i  
  }//k

}

#endif /* TRACKS_H */
