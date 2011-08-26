// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef TRACKS_H
#define TRACKS_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

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

namespace vw
{
namespace math
{
  class pointCloud : public Vector<double,3>
  {
  public:
  
    pointCloud
      (
  	Vector3 = Vector3(), //coords
        int   = 0, // year
  	int   = 0, // month
  	int   = 0, // day
  	int   = 0, // hour
  	int   = 0, // min
  	float = 0, // sec
  	int   = 0  // s
      );
  
    ~pointCloud(){}; 
  
    int year;
    int month;
    int day;
    int hour; //24hr format
    int min;
    float sec;
    int s; //detector id
  };
  
  inline
  std::ostream& operator<< ( std::ostream& stream, pointCloud p )
    {
    stream 
      << p.year
      << "-"
      << p.month
      << "-"
      << p.day
      << "T"
      << p.hour
      << ":"
      << p.min
      << ":"
      << p.sec
      << " S: "
      << p.s
      << "  Vector( " << p.x() << ", " << p.y() << ", " << p.z() << " )" ;
  
    return stream;
    }
  
}
}

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

class LOLAShot
{
  public:

  explicit LOLAShot( vector<pointCloud> );
  explicit LOLAShot( pointCloud );

  ~LOLAShot(){};

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

  private:
  void init
    (
    vector<pointCloud>,
    vector<imgPoint> = vector<imgPoint>(), // explicitly make empty
    vector<DEMPoint> = vector<DEMPoint>(), // explicitly make empty
    int              = 0, // valid
    int              = 0, // centerPtIndex
    float            = 0, // reflectance
    float            = 0, // synthImage
    int              = 0, // calc_acp
    float            = 0, // filter_response
    float            = 0, // featurePtRefl
    float            = 1, // weightRefl
    float            = 0, // featurePtLOLA
    float            = 0  // filresLOLA
    );
};

inline
std::ostream& operator<< ( std::ostream& stream, LOLAShot s )
  {
  stream 
    << "valid: " << s.valid
    << " centerPtIndex: " << s.centerPtIndex
    << " reflectance: "   << s.reflectance
    << " synthiImage: "   << s.synthImage
    << " calc_acp: "      << s.calc_acp
    << endl
    << "featurePtRefl: "  << s.featurePtRefl
    << " weightRefl: "    << s.weightRefl
    << " featurePtLOLA: " << s.featurePtLOLA
    << " filresLOLA: "    << s.filresLOLA
    << endl
    << "Number of LOLA points in this shot: " << s.LOLAPt.size();
  for( unsigned int point = 0; point < s.LOLAPt.size(); ++point )
    { stream << endl << "  " << s.LOLAPt[point]; }

  stream << endl << "Number of img points in this shot: " << s.imgPt.size();
  for( unsigned int point = 0; point < s.imgPt.size(); ++point )
    {
    stream << endl
           << " x: "   << s.imgPt[point].x 
           << " y: "   << s.imgPt[point].y 
           << " val: " << s.imgPt[point].val;
    }

  stream << endl << "Number of DEM points in this shot: " << s.DEMPt.size();
  for( unsigned int point = 0; point < s.DEMPt.size(); ++point )
    {
    stream << endl
           << " valid: " << s.DEMPt[point].valid 
           << " x: "     << s.DEMPt[point].x
           << " y: "     << s.DEMPt[point].y
           << " val: "   << s.DEMPt[point].val;
    }

  return stream;
  }


int GetTimeDiff(pointCloud prevPt, pointCloud currPt, float timeThresh);
vector<vector<LOLAShot> > CSVFileRead(string CSVFilename);
//vector<vector<LOLAShot> > CSVFileRead_LIMA(string CSVFilename);
//computes the scale factor for all tracks at once
Vector4 FindMinMaxLat(vector<vector<LOLAShot> >trackPts);

//float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts);
Vector2 ComputeGainBiasFactor(vector<vector<LOLAShot > >&trackPts);
int GetAllPtsFromCub(vector<vector<LOLAShot > > &trackPts, string cubFilename);
int ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks, ModelParams modelParams, CoregistrationParams coregistrationParams);
int ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks,  Vector3 cameraPosition, Vector3 lightPosition, CoregistrationParams coregistrationParams);
pointCloud GetPointFromIndex(vector<pointCloud> const &  LOLAPts, int index);

void SaveReflectancePoints(vector< vector<LOLAShot> >  &allTracks, Vector2 gain_bias, string filename);
void SaveImagePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename);
void SaveAltitudePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename);
void SaveDEMPoints(vector< vector<LOLAShot> > &tracks, string DEMFilename, string filename);

void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray,  string cubFile, 
               vector<gcp> &gcpArray, Vector2 centroid);

void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray, 
               string camCubFile, string mapCubFile, vector<gcp> &gcpArray, Vector2 centroid, 
               float downsample_factor);
void UpdateGCP(vector<vector<LOLAShot> > trackPts, 
               vector<Vector4> matchArray, vector<float> errorArray,
               string camCubFile, string mapCubFile, 
               vector<gcp> &gcpArray, float downsample_factor);

void SaveGCPoints(vector<gcp> gcpArray,  string gcpFilename);

vector<float> GetTrackPtsByID(vector<LOLAShot> trackPts, int ID);
vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID);
void ComputeAverageShotDistance(vector<vector<LOLAShot> >trackPts);
void ComputeAverageIntraShotDistance(vector<vector<LOLAShot> >trackPts);
//this function does not belong here
Vector2 ComputeMinMaxValuesFromCub(string cubFilename);



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
  } 
  else 
  {
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
          
	    float lon = LOLAPts[j].x();
	    float lat = LOLAPts[j].y();
	    float rad = LOLAPts[j].z();
     
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

//determines LOLA corresponding points in a DEM
template <class ViewT>
void 
GetAllPtsFromDEM(vector<vector<LOLAShot> >&  trackPts,
                 const ImageViewBase<ViewT>& DEM,
                 const GeoReference&         DEMGeo,
                 const double                noDEMVal)
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


  for(unsigned int ti = 0; ti < trackPts.size(); ti++){
    for(unsigned int si = 0; si < trackPts[ti].size(); si++){
   
      LOLAPts = trackPts[ti][si].LOLAPt;
  
      trackPts[ti][si].valid = 0;
      if (LOLAPts.size() > 0){ 
        trackPts[ti][si].valid = 1;
	trackPts[ti][si].DEMPt.resize(LOLAPts.size());
	
	for (unsigned int li = 0; li < LOLAPts.size(); li++){
	  float lon = LOLAPts[li].x();
	  float lat = LOLAPts[li].y();
	  // float rad = LOLAPts[li].coords[2];
	  
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
      } 
    }//i  
  }//k
}


#endif /* TRACKS_H */
