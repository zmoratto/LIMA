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
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
#include "coregister.h"

struct ShotTime
{
  int year;
  int month;
  int dap;
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

struct DEMPoint
{
  int valid;
  double x;
  double y;
  double val;
};


// prototype for GetPointFromIndex
pointCloud GetPointFromIndex( vector<pointCloud> const & LOLAPts, int index);

struct LOLAShot
{
  int valid;
  int centerPtIndex;
  float reflectance;
  float synthImage;

  //following variables are for LOLA feature and weight computation
  int   calc_acp;             //is the filter valid here?
  float filter_response;     
  float featurePtRefl;          
  float weightRefl;        


  float featurePtLOLA;
  float filresLOLA;           
        
  vector<pointCloud> LOLAPt;
  vector<imgPoint> imgPt;
  vector<DEMPoint> DEMPt; 
};

int GetTimeDiff(pointCloud prevPt, pointCloud currPt, float timeThresh);
//computes the scale factor for all tracks at once
Vector4 FindMinMaxLat(vector<vector<LOLAShot> >trackPts);
vector<vector<LOLAShot> > CSVFileRead(string CSVFilename);
Vector2 ComputeMinMaxValuesFromCub(string cubFilename);
float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts);
void  ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks, ModelParams modelParams, CoregistrationParams coregistrationParams);
pointCloud GetPointFromIndex(vector<pointCloud> const &  LOLAPts, int index);


void SaveReflectancePoints(vector< vector<LOLAShot> >  &allTracks, float scaleFactor, string filename);
void SaveImagePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename);
void SaveAltitudePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename);
void SaveDEMPoints(vector< vector<LOLAShot> > &tracks, string DEMFilename, string filename);

void SaveGCPoints(vector<vector<LOLAShot> > trackPts,  std::vector<std::string> DRGFiles,  std::vector<int> overlapIndices, 
                  vector<Vector<float, 6> > optimalTransfArray, vector<float> optimalErrorArray, string gcpFilename);

vector<float> GetTrackPtsByID(vector<LOLAShot> trackPts, int ID);
vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID);

//template <class ViewT>
void 
GetAllPtsFromCub(vector<vector<LOLAShot > > &trackPts, string cubFilename);

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

    cout << "NOT masked" <<endl;
  } else {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(apply_mask(DRG.impl()),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );
    cout << "MASKED" <<endl;
  }

  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      
      trackPts[k][i].valid = 1; 

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

            //check that (x,y) are within the image boundaries
	    if ((x>=0) && (y>=0) && (x<DRG.impl().cols()) && (y<DRG.impl().rows())){//valid position  
              //check for valid data as well
              if (interpDRG(x, y)!=0){//valid values
	         trackPts[k][i].imgPt[j].val = interpDRG(x, y);
	         trackPts[k][i].imgPt[j].x = DRG_pix[0];
	         trackPts[k][i].imgPt[j].y = DRG_pix[1];
	      }
              else{//invalidate the point
                 trackPts[k][i].valid = 0;
              }
	    }
            else{ //invalidate the point  
                 trackPts[k][i].valid = 0; 
            }
	
      }
      //check for valid shot with 3 points to compute reflectance
      pointCloud centerPt  = GetPointFromIndex( LOLAPts, 3);
      pointCloud topPt     = GetPointFromIndex( LOLAPts, 2);
      pointCloud leftPt    = GetPointFromIndex( LOLAPts, 1);
      if ((centerPt.s == -1) || (topPt.s == -1) || (leftPt.s == -1) || (LOLAPts.size() >5)){//invalid LOLA shot
          trackPts[k][i].valid = 0; 
      }

    }//i  
  }//k

}


template <class ViewT>
void 
GetAllPtsFromDEM(vector<vector<LOLAShot > > &trackPts,  ImageViewBase<ViewT> const& DEM, GeoReference const &DEMGeo, double noDEMVal)
{
  vector<pointCloud> LOLAPts;

  float radius = DEMGeo.datum().semi_major_axis();
 
  ImageViewRef<float> interpDEM;
  if ( IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDEM = pixel_cast<float>(interpolate(edge_extend(DEM.impl(),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );

    //cout << "NOT masked" <<endl;
  } else {
    interpDEM = pixel_cast<float>(interpolate(edge_extend(apply_mask(DEM.impl()),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }


  for(int ti = 0; ti < trackPts.size(); ti++){
    for(int si = 0; si < trackPts[ti].size(); si++){
   
      LOLAPts = trackPts[ti][si].LOLAPt;
      int numLOLAPts = LOLAPts.size();
  
      trackPts[ti][si].valid = 0;
      if (numLOLAPts > 0){ 

	trackPts[ti][si].DEMPt.resize(LOLAPts.size());

	for (int li = 0; li < LOLAPts.size(); li++){
          
	    float lon = LOLAPts[li].coords[0];
	    float lat = LOLAPts[li].coords[1];
	    float rad = LOLAPts[li].coords[2];
     
	    Vector2 DEM_lonlat(lon, lat);
	    Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(DEM_lonlat);
	  
	    float x = DEM_pix[0];
	    float y = DEM_pix[1];
      
            trackPts[ti][si].DEMPt[li].valid = 0; 
            if ((x>=0) && (y>=0) && (x<interpDEM.cols()) && (y<interpDEM.rows())){
	      if (interpDEM(x,y)!=noDEMVal){
		trackPts[ti][si].DEMPt[li].val = 0.001*(radius + interpDEM(x, y));
		trackPts[ti][si].DEMPt[li].x = DEM_pix[0];
		trackPts[ti][si].DEMPt[li].y = DEM_pix[1];
                trackPts[ti][si].DEMPt[li].valid=1; 
	      }
	    }
	   
	}

	trackPts[ti][si].valid = 1;
      } 
     
    }//i  
  }//k
}



#endif /* TRACKS_H */
