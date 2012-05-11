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
#include <vw/Math.h>
#include <vw/Math/BBox.h>
#include <vw/Math/Vector.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;
#include <math.h>
#include "coregister.h"

class pointCloud : public vw::Vector<double,3> {
  public:
  
  pointCloud( const vw::Vector3& = vw::Vector3(), //coords
              const int&         = 0,             // year
              const int&         = 0,             // month
              const int&         = 0,             // day
              const int&         = 0,             // hour
              const int&         = 0,             // min
              const float&       = 0,             // sec
              const int&         = 0              // s
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
pointCloud GetPointFromIndex( const std::vector<pointCloud>&, const int );

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


std::vector<std::vector<LOLAShot> > LOLAFileRead( const std::string& );
inline std::vector<std::vector<LOLAShot> > CSVFileRead( const std::string& f ){
  // This function is deprecated, please use LOLAFileRead().
  return LOLAFileRead( f );
}

// Finds the lat/lon bounding box that encloses all the points.
vw::BBox2  FindShotBounds( const std::vector<std::vector<LOLAShot> >& );
vw::Vector4 FindMinMaxLat( const std::vector<std::vector<LOLAShot> >& );


//float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts);
vw::Vector2 ComputeGainBiasFactor( const std::vector<LOLAShot>& );
vw::Vector2 ComputeGainBiasFactor( const std::vector<std::vector<LOLAShot> >& );
int GetAllPtsFromCub(vector<vector<LOLAShot > > &trackPts, string cubFilename);
int ComputeAllReflectance(       std::vector< std::vector<LOLAShot> >& shots,
                           const vw::Vector3&                          cameraPosition,
                           const vw::Vector3&                          lightPosition);


void SaveReflectancePoints( const std::vector<std::vector<LOLAShot> >&,
                            const vw::Vector2&                        gain_bias,
                            const std::string&                        filename);

void SaveImagePoints( const std::vector<std::vector<LOLAShot> >&,
                      const int&                                detectNum, 
                      const std::string&                        filename);

void SaveAltitudePoints( const std::vector<std::vector<LOLAShot> >&,
                         const int&                                detectNum,
                         const std::string&                        filename);


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

// These are not currently being used.
// void ComputeAverageShotDistance(vector<vector<LOLAShot> >trackPts);
// void ComputeAverageIntraShotDistance(vector<vector<LOLAShot> >trackPts);

//this function does not belong here
vw::Vector2 ComputeMinMaxValuesFromCub( const std::string& );
//Vector2 ComputeMinMaxValuesFromDEM(string demFilename);


template <class ViewT>
void 
GetAllPtsFromImage(       std::vector<std::vector<LOLAShot> >& trackPts,  
                    const vw::ImageViewBase<ViewT>&            DRG, 
                    const GeoReference&                        DRGGeo) {
  ImageViewRef<float> interpDRG;
  if( IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(DRG.impl(),
                                                          ConstantEdgeExtension()),
                                              BilinearInterpolation()) );
    //cout << "NOT masked" <<endl;
  } 
  else {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(apply_mask(DRG.impl()),
                                                          ConstantEdgeExtension()),
                                              BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }

  for( unsigned int k = 0; k < trackPts.size(); ++k ){
    for( unsigned int i = 0; i < trackPts[k].size(); ++i ){

      // Set initially valid
      trackPts[k][i].valid = 1; 

      trackPts[k][i].imgPt.resize( trackPts[k][i].LOLAPt.size() );

      BBox2i bbox = bounding_box( DRG );
      for( unsigned int j = 0; j < trackPts[k][i].LOLAPt.size(); ++j ){
	    Vector2 DEM_lonlat( trackPts[k][i].LOLAPt[j].x(), trackPts[k][i].LOLAPt[j].y() );
	    Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);
	  
	    if( bbox.contains(DRG_pix) ){
          //check for valid data as well
          if( interpDRG( DRG_pix[0], DRG_pix[1] ) != 0 ){//valid values
            trackPts[k][i].imgPt[j].val = interpDRG( DRG_pix[0], DRG_pix[1] );
            trackPts[k][i].imgPt[j].x = DRG_pix[0];
            trackPts[k][i].imgPt[j].y = DRG_pix[1];
          }
          else{// one bad point invalidates the whole set of shots
            trackPts[k][i].valid = 0;
          }
        }
        else { trackPts[k][i].valid = 0; }
      }
      //check for valid shot with 3 points to compute reflectance
      pointCloud centerPt  = GetPointFromIndex( trackPts[k][i].LOLAPt, 3);
      pointCloud topPt     = GetPointFromIndex( trackPts[k][i].LOLAPt, 2);
      pointCloud leftPt    = GetPointFromIndex( trackPts[k][i].LOLAPt, 1);
      if ((centerPt.s == -1) || (topPt.s == -1) || (leftPt.s == -1) || (trackPts[k][i].LOLAPt.size() >5)){//invalid LOLA shot
          trackPts[k][i].valid = 0; 
      }
    }//i  
  }//k

}

//determines LOLA corresponding points in a DEM
template <class ViewT>
void GetAllPtsFromDEM(       std::vector<std::vector<LOLAShot> >& trackPts,
                       const vw::ImageViewBase<ViewT>&            DEM,
                       const vw::cartography::GeoReference&       DEMGeo,
                       const double&                              noDEMVal) {
  //determine the minmx value of the DEM - START
  float minVal = std::numeric_limits<float>::max();
  float maxVal = std::numeric_limits<float>::min();
  min_max_pixel_values( create_mask(DEM,noDEMVal), minVal, maxVal );
  //cout<<"min="<<minVal<<", max="<<maxVal<<endl; 
  //determine the minmx value of the DEM - END

  vw::ImageViewRef<float> interpDEM;
  if( vw::IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDEM = vw::pixel_cast<float>(interpolate(edge_extend(DEM.impl(),
							          vw::ConstantEdgeExtension()),
					                  vw::BilinearInterpolation()) );
    //cout << "NOT masked" <<endl;
  } else {
  interpDEM = vw::pixel_cast<float>(interpolate(edge_extend(apply_mask(DEM.impl()),
							      vw::ConstantEdgeExtension()),
					              vw::BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }

  vw::BBox2i iDEM_bbox = vw::bounding_box( interpDEM );

  float radius = DEMGeo.datum().semi_major_axis();
  for( unsigned int ti = 0; ti < trackPts.size(); ++ti ){
    for( unsigned int si = 0; si < trackPts[ti].size(); ++si ){
//      LOLAPts = trackPts[ti][si].LOLAPt;
      trackPts[ti][si].valid = 0;
      if( trackPts[ti][si].LOLAPt.size() > 0 ){ 
        trackPts[ti][si].valid = 1;
        trackPts[ti][si].DEMPt.resize( trackPts[ti][si].LOLAPt.size() );
	
        for( unsigned int li = 0; li < trackPts[ti][si].LOLAPt.size(); ++li ){
          vw::Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(
                                       subvector( trackPts[ti][si].LOLAPt[li], 0, 2 ) );
          float x = DEM_pix.x();
          float y = DEM_pix.y();
      
          trackPts[ti][si].DEMPt[li].valid = 0; 
          if( (iDEM_bbox.contains(DEM_pix)) && 
              (interpDEM(x,y) > minVal)     && 
              (interpDEM(x,y) < maxVal)       ){
            trackPts[ti][si].DEMPt[li].x = x;
            trackPts[ti][si].DEMPt[li].y = y;
            trackPts[ti][si].DEMPt[li].val = 0.001*(radius + interpDEM(x, y));
            trackPts[ti][si].DEMPt[li].valid = 1; 
          }
        }
      }
    } 
  }
}

//determines LOLA corresponding points in a DEM
template <class ViewT>
void GetAllPtsFromDEM_Prec(       std::vector<std::vector<LOLAShot> >& trackPts,
                            const vw::ImageViewBase<ViewT>&            DEM,
                            const vw::cartography::GeoReference&       DEMGeo,
                            const double&                              noDEMVal,
                            const vw::ImageViewBase<ViewT>&            DEM_Prec) {
  //determine the minmx value of the DEM - START
  float minVal = std::numeric_limits<float>::max();
  float maxVal = std::numeric_limits<float>::min();
  min_max_pixel_values( create_mask(DEM,noDEMVal), minVal, maxVal );
  //cout<<"min="<<minVal<<", max="<<maxVal<<endl; 
  //determine the minmx value of the DEM - END

  vw::ImageViewRef<float> interpDEM;
  if( vw::IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDEM = vw::pixel_cast<float>(interpolate(edge_extend(DEM.impl(),
                                      vw::ConstantEdgeExtension()),
                                      vw::BilinearInterpolation()) );
    //cout << "NOT masked" <<endl;
  } else {
    interpDEM = vw::pixel_cast<float>(interpolate(edge_extend(apply_mask(DEM.impl()),
                                      vw::ConstantEdgeExtension()),
                                      vw::BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }

  vw::ImageViewRef<float> interpDEM_Prec;
  interpDEM_Prec = vw::pixel_cast<float>(interpolate(edge_extend(DEM_Prec.impl(),
                                         vw::ConstantEdgeExtension()),
                                         vw::BilinearInterpolation()) );

  vw::BBox2i iDEM_bbox = vw::bounding_box( interpDEM );

  float radius = DEMGeo.datum().semi_major_axis();
  for( unsigned int ti = 0; ti < trackPts.size(); ++ti ){
    for( unsigned int si = 0; si < trackPts[ti].size(); ++si ){
      std::vector<pointCloud> points = trackPts[ti][si].LOLAPt;

      trackPts[ti][si].valid = 0;
      if( points.size() > 0 ){ 
        trackPts[ti][si].valid = 1;
        trackPts[ti][si].DEMPt.resize( points.size() );
	
        for( unsigned int li = 0; li < points.size(); ++li ){
          vw::Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(
                                       subvector( points[li], 0, 2 ) );
          float x = DEM_pix.x();
          float y = DEM_pix.y();
      
          trackPts[ti][si].DEMPt[li].valid = 0; 
          if( (iDEM_bbox.contains(DEM_pix)) &&
              (interpDEM(x,y) > minVal)     && 
              (interpDEM(x,y) < maxVal)     && 
              (interpDEM_Prec(x,y) >= 0)    && 
              (interpDEM_Prec(x,y) < 25)      ){

            trackPts[ti][si].DEMPt[li].x = x;
            trackPts[ti][si].DEMPt[li].y = y;
            trackPts[ti][si].DEMPt[li].val = 0.001*(radius + interpDEM(x, y));
            trackPts[ti][si].DEMPt[li].valid = 1; 
          }
        }
      }
    }
  }
}

#endif /* TRACKS_H */

