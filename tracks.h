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

#include <asp/IsisIO/IsisCameraModel.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;
#include <math.h>
#include "coregister.h"

// from a lidar track, 3D points and time of capture
class pointCloud : public vw::Vector<double,3> {
  public:
  
  pointCloud( const vw::Vector3& = vw::Vector3(), //coords
              const int          = 0,             // year
              const int          = 0,             // month
              const int          = 0,             // day
              const int          = 0,             // hour
              const int          = 0,             // min
              const float        = 0,             // sec
              const int          = 0              // s
              );
  
  ~pointCloud(){};

  inline std::string time(){
    ostringstream os;
    os << this->year
       << "-"
       << setw(2) << setfill('0') << this->month
       << "-"
       << setw(2) << setfill('0') << this->day
       << "T"
       << setw(2) << setfill('0') << this->hour
       << ":"
       << setw(2) << setfill('0') << this->min
       << ":";
    if( this->sec < 10 ){ os << "0"; }
    os << this->sec;
    return os.str();
  }
  
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
      << p.time()
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

// regular points on a surface mesh
struct DEMPoint
{
  int valid;// map may have holes in the terrain
  double x;
  double y;
  double val;
};

bool isTimeDiff( const pointCloud&, const pointCloud&, const float);

// vector of point clouds, 5 points (but sometimes more or less)
// prototype for GetPointFromIndex
pointCloud GetPointFromIndex( const std::vector<pointCloud>&, const int );

class LOLAShot {
  public:

  explicit LOLAShot( const std::vector<pointCloud>& );
  explicit LOLAShot( const pointCloud& );

  ~LOLAShot(){};

  int valid;
  // int centerPtIndex;  Not used?
  float reflectance;

  //following variables are for LOLA feature and weight computation
  int   calc_acp;             //is the filter valid here?
  float filter_response;     
  float featurePtRefl;          
  float weightRefl;        


  float featurePtLOLA;
  float filresLOLA;

  float image;// the  color in the original image
        
  std::vector<pointCloud> LOLAPt;
  // points on image and DEM corresponding to point cloud
  std::vector<imgPoint>   imgPt; // original locations of each shot in the original image
  std::vector<DEMPoint>   DEMPt; 

  private:
  void init
    (
    std::vector<pointCloud>,
    std::vector<imgPoint> = std::vector<imgPoint>(), // explicitly make empty
    std::vector<DEMPoint> = std::vector<DEMPoint>(), // explicitly make empty
    int              = 0, // valid
    // int              = 0, // centerPtIndex
    float            = 0, // reflectance
    // float            = 0, // synthImage
    int              = 0, // calc_acp
    float            = 0, // filter_response
    float            = 0, // featurePtRefl
    float            = 1, // weightRefl
    float            = 0, // featurePtLOLA
    float            = 0  // filresLOLA
    );
};

class AlignedLOLAShot : public LOLAShot
{
public:
  explicit AlignedLOLAShot( LOLAShot& s) : LOLAShot(s), image_x(-1), image_y(-1), synth_image(-1) {};
  ~AlignedLOLAShot(){};

  int image_x, image_y; // adjusted position in the image
  float synth_image; // predicted color after applying gain

};

inline
std::ostream& operator<< ( std::ostream& stream, LOLAShot s )
  {
  stream 
    << "valid: " << s.valid
    // << " centerPtIndex: " << s.centerPtIndex
    << " reflectance: "   << s.reflectance
    // << " synthiImage: "   << s.synthImage
    << " calc_acp: "      << s.calc_acp
    << std::endl
    << "featurePtRefl: "  << s.featurePtRefl
    << " weightRefl: "    << s.weightRefl
    << " featurePtLOLA: " << s.featurePtLOLA
    << " filresLOLA: "    << s.filresLOLA
    << std::endl
    << "Number of LOLA points in this shot: " << s.LOLAPt.size();
  for( unsigned int point = 0; point < s.LOLAPt.size(); ++point )
    { stream << std::endl << "  " << s.LOLAPt[point]; }

  stream << endl << "Number of img points in this shot: " << s.imgPt.size();
  for( unsigned int point = 0; point < s.imgPt.size(); ++point )
    {
    stream << std::endl
           << " x: "   << s.imgPt[point].x 
           << " y: "   << s.imgPt[point].y 
           << " val: " << s.imgPt[point].val;
    }

  stream << std::endl << "Number of DEM points in this shot: " << s.DEMPt.size();
  for( unsigned int point = 0; point < s.DEMPt.size(); ++point )
    {
    stream << std::endl
           << " valid: " << s.DEMPt[point].valid 
           << " x: "     << s.DEMPt[point].x
           << " y: "     << s.DEMPt[point].y
           << " val: "   << s.DEMPt[point].val;
    }

  return stream;
  }


std::vector<std::vector<LOLAShot> > LOLAFileRead( const std::string& );
std::vector<std::vector<LOLAShot> > LOLAFilesRead( const std::string& );
inline std::vector<std::vector<LOLAShot> > CSVFileRead( const std::string& f ){
  // This function is deprecated, please use LOLAFileRead().
  return LOLAFileRead( f );
}

// Finds the lat/lon bounding box that encloses all the points.
vw::BBox2  FindShotBounds( const std::vector<std::vector<LOLAShot> >& );
vw::Vector4 FindMinMaxLat( const std::vector<std::vector<LOLAShot> >& );


//float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts);
vw::Vector2 ComputeGainBiasFactor( const std::vector<AlignedLOLAShot>& );
vw::Vector2 ComputeGainBiasFactor( const std::vector<std::vector<AlignedLOLAShot> >& );
int GetAllPtsFromCub( vector<vector<LOLAShot> >& trackPts, camera::IsisCameraModel & model, ImageView<PixelGray<float> > & cubImage);
int ComputeAllReflectance(       std::vector< std::vector<LOLAShot> >& shots,
                           const vw::Vector3&                          cameraPosition,
                           const vw::Vector3&                          lightPosition);


void SaveReflectancePoints( const std::vector<std::vector<LOLAShot> >&,
                            const vw::Vector2&                        gain_bias,
                            const std::string&                        filename);

void SaveImagePoints( const std::vector<std::vector<LOLAShot> >&,
                      const int                                 detectNum, 
                      const std::string&                        filename);

void SaveAltitudePoints( const std::vector<std::vector<LOLAShot> >&,
                         const int                                 detectNum,
                         const std::string&                        filename);

/* These versions aren't used.  Are they deprecated?
* void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray,  string cubFile, 
*                vector<gcp> &gcpArray, Vector2 centroid);
* 
* void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray, 
*                string camCubFile, string mapCubFile, vector<gcp> &gcpArray, Vector2 centroid, 
*                float downsample_factor);
*/
void UpdateGCP( const std::vector<std::vector<LOLAShot> >& trackPts, 
                const std::vector<vw::Vector4>& matchArray, 
                const std::string& camCubFile, 
                const std::string& mapCubFile, 
                std::vector<gcp>& gcps,
                float downsample_factor);

void UpdateGCP( const vector<vector<AlignedLOLAShot> >& trackPts, 
                const string&                    camCubFile, 
                vector<gcp>&                     gcpArray);

void SaveGCPoints( const std::vector<gcp>&, const std::string& );

// These are not currently being used.
// void ComputeAverageShotDistance(vector<vector<LOLAShot> >trackPts);
// void ComputeAverageIntraShotDistance(vector<vector<LOLAShot> >trackPts);

//this function does not belong here
vw::Vector2 ComputeMinMaxValuesFromCub( ImageView<PixelGray<float> > &  image );
//Vector2 ComputeMinMaxValuesFromDEM(string demFilename);


template <class ViewT>
void 
GetAllPtsFromImage(       std::vector<std::vector<LOLAShot> >& trackPts,  
                    const vw::ImageViewBase<ViewT>&            DRG, 
                    const vw::cartography::GeoReference&       DRGGeo) {
  ImageViewRef<float> interpDRG;
  if( IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDRG = vw::pixel_cast<float>(interpolate(edge_extend(DRG.impl(),
                                                          vw::ConstantEdgeExtension()),
                                              vw::BilinearInterpolation()) );
    //cout << "NOT masked" <<endl;
  } 
  else {
    interpDRG = vw::pixel_cast<float>(interpolate(edge_extend(apply_mask(DRG.impl()),
                                                          vw::ConstantEdgeExtension()),
                                              vw::BilinearInterpolation()) );
    //cout << "MASKED" <<endl;
  }

  for( unsigned int k = 0; k < trackPts.size(); ++k ){
    for( unsigned int i = 0; i < trackPts[k].size(); ++i ){
      vector<pointCloud> points = trackPts[k][i].LOLAPt;
      if( points.size() >5){ 
        trackPts[k][i].valid = 0;
        break;
      }

      // Set initially valid
      trackPts[k][i].valid = 1; 

      trackPts[k][i].imgPt.resize( points.size() );

      BBox2i bbox = bounding_box( DRG );
      for( unsigned int j = 0; j < points.size(); ++j ){
	    Vector2 DEM_lonlat( points[j].x(), points[j].y() );
	    Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);
	  
	    if( (bbox.contains(DRG_pix))
            && (interpDRG( DRG_pix.x(), DRG_pix.y() ) != 0 )){
          trackPts[k][i].imgPt[j].val = interpDRG( DRG_pix.x(), DRG_pix.y() );
          trackPts[k][i].imgPt[j].x = DRG_pix.x();
          trackPts[k][i].imgPt[j].y = DRG_pix.y();
        }
        else{// one bad point invalidates the whole set of shots
          trackPts[k][i].valid = 0;
          break;
        }
      }
      //check for valid shot with 3 points to compute reflectance
      try {
        GetPointFromIndex( points, 3);
        GetPointFromIndex( points, 2);
        GetPointFromIndex( points, 1);
      }
      catch( const vw::ArgumentErr& error ){ trackPts[k][i].valid = 0;}
    }//i  
  }//k

}

//determines LOLA corresponding points in a DEM
template <class ViewT>
void GetAllPtsFromDEM(       std::vector<std::vector<LOLAShot> >& trackPts,
                       const vw::ImageViewBase<ViewT>&            DEM,
                       const vw::cartography::GeoReference&       DEMGeo,
                       const double                               noDEMVal) {
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
      vector<pointCloud> points = trackPts[ti][si].LOLAPt;
      trackPts[ti][si].valid = 0;
      if( points.size() > 0 ){ 
        trackPts[ti][si].valid = 1;
        trackPts[ti][si].DEMPt.resize( points.size() );
	
        for( unsigned int li = 0; li < points.size(); ++li ){
          vw::Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(
                                       subvector( points[li], 0, 2 ) );
          trackPts[ti][si].DEMPt[li].valid = 0; 
          if( (iDEM_bbox.contains(DEM_pix))                    && 
              (interpDEM( DEM_pix.x(), DEM_pix.y() ) > minVal) && 
              (interpDEM( DEM_pix.x(), DEM_pix.y() ) < maxVal) ){
            trackPts[ti][si].DEMPt[li].x = DEM_pix.x();
            trackPts[ti][si].DEMPt[li].y = DEM_pix.y();
            trackPts[ti][si].DEMPt[li].val = 0.001*(radius + interpDEM( DEM_pix.x(), DEM_pix.y() ));
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
                            const double                               noDEMVal,
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

void transform_track(vector<AlignedLOLAShot> & track, Matrix3x3 transform, ImageView<PixelGray<float> >& cub);
void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Matrix3x3 transform, ImageView<PixelGray<float> >& cub);
void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Matrix3x3 transform, string cubFile);
void transform_tracks_by_matrices(vector<vector<LOLAShot> > & tracks, vector<Matrix3x3> matrices);
void transform_tracks_by_matrix(vector<vector<LOLAShot> > & tracks, Matrix3x3 M);

void save_track_data( const std::vector<std::vector<AlignedLOLAShot> >&, const std::string& filename);
vector<Matrix3x3> load_track_transforms(const std::string& filename);

std::vector<std::vector< AlignedLOLAShot> > initialize_aligned_lola_shots(std::vector<std::vector<LOLAShot> >& tracks);

#endif /* TRACKS_H */

