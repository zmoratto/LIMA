#include <boost/filesystem.hpp>

#include "../tracks.h"
#include "../icp.h"
#include "gtest/gtest.h"

namespace fs = boost::filesystem;

using namespace std;

TEST( ComputeDEMTranslation, works ){
  vector<Vector3> features;
  features.push_back( Vector3(1,1,1) );
  features.push_back( Vector3(2,3,4) );
  features.push_back( Vector3(21,31,41) );

  vector<Vector3> match;
  match.push_back( Vector3(0,0,0) );
  match.push_back( Vector3(1,2,3) );
  match.push_back( Vector3(20,30,40) );

  Vector3 truth( -1, -1, -1 );

  Vector3 translated = ComputeDEMTranslation( features, match );
  ASSERT_EQ( truth, translated ) << "The translation is wrong.";
}

class ComputeMatchingError_Test : public ::testing::Test {
  protected:
  virtual void SetUp(){
    features_.push_back( Vector3(1,1,1) );
    features_.push_back( Vector3(2,3,4) );
    features_.push_back( Vector3(21,31,41) );
    match_.push_back( Vector3(0,0,0) );
    match_.push_back( Vector3(1,2,3) );
    match_.push_back( Vector3(20,30,40) );
  }

  virtual void TearDown() {}

  vector<Vector3> features_;
  vector<Vector3> match_;
};

TEST_F( ComputeMatchingError_Test, basic ){
  valarray<float> truth( sqrt(3), 3 );
  valarray<float> errors = ComputeMatchingError( features_, match_ );

  for( unsigned int i = 0; i < errors.size(); ++i ){
    EXPECT_NEAR( truth[i], errors[i], 0.001 ) << "Error array mismatch.";
  }
}

TEST_F( ComputeMatchingError_Test, 3D){
  float error = ComputeMatchingError3D( features_, match_ );
  ASSERT_NEAR( sqrt(3), error, 0.001 ) << "Incorrect error.";
}


class ComputeDEMRotation_Test : public ::testing::Test {
  protected:
  virtual void SetUp(){
    features_.push_back( Vector3( 5,5,6 ) );
    features_.push_back( Vector3( 5,6,5 ) );
    features_.push_back( Vector3( 6,5,5 ) );

    match_.push_back( Vector3( 5,5,4 ) );
    match_.push_back( Vector3( 5,6,5 ) );
    match_.push_back( Vector3( 6,5,5 ) );

    rotation_(0,0) =  1.0/3.0;
    rotation_(0,1) = -2.0/3.0;
    rotation_(0,2) = -2.0/3.0;
    rotation_(1,0) = -2.0/3.0;
    rotation_(1,1) =  1.0/3.0;
    rotation_(1,2) = -2.0/3.0;
    rotation_(2,0) =  2.0/3.0;
    rotation_(2,1) =  2.0/3.0;
    rotation_(2,2) = -1.0/3.0;
  }

  virtual void TearDown() {}
  
  vector<Vector3> features_;
  vector<Vector3> match_;
  Matrix<float, 3, 3> rotation_;
};
  
TEST_F( ComputeDEMRotation_Test, 4arg ){

  Vector3 f_center = find_centroid( features_ );
  Vector3 m_center = find_centroid( match_ );
  Matrix<float, 3, 3> test = ComputeDEMRotation( features_, match_, f_center, m_center );

  for( unsigned int i = 0; i < 3; ++i ){
    for( unsigned int j = 0; j < 3; ++j ){
    EXPECT_NEAR( rotation_(i,j), test(i,j), 0.0001 ) << "Rotation elements aren't correct.";
    }
  }
}

TEST_F( ComputeDEMRotation_Test, 2arg ){

  Matrix<float, 3, 3> test = ComputeDEMRotation( features_, match_ );

  for( unsigned int i = 0; i < 3; ++i ){
    for( unsigned int j = 0; j < 3; ++j ){
    EXPECT_NEAR( rotation_(i,j), test(i,j), 0.0001 ) << "Rotation elements aren't correct.";
    }
  }
}

TEST_F( ComputeDEMRotation_Test, TransformFeatures4arg ){
  
  Vector3 translation( 0,0,0 );
  Vector3 f_Center = find_centroid( features_ );
  TransformFeatures( features_, f_Center, translation, rotation_ );

  for( unsigned int i = 0; i < features_.size(); ++i ){
    Vector3 truth = match_[i] + Vector3( 0,0,0.6666 ) + translation;
    for( unsigned int j = 0; j < 3; ++j ){
      EXPECT_NEAR( truth(j), features_[i](j), 0.0001 ) << "Translated and rotated elements aren't correct.";
    }
  }
}

TEST_F( ComputeDEMRotation_Test, TransformFeatures3arg ){
  
  Vector3 translation( 1,2,3 );
  TransformFeatures( features_, translation, rotation_ );

  for( unsigned int i = 0; i < features_.size(); ++i ){
    Vector3 truth = match_[i] + Vector3( 0,0,0.6666 ) + translation;
    for( unsigned int j = 0; j < 3; ++j ){
      EXPECT_NEAR( truth(j), features_[i](j), 0.0001 ) << "Translated and rotated elements aren't correct.";
    }
  }
}

class Features_Test : public ::testing::Test {
  protected:
  virtual void SetUp(){
  fives_.x() = 5;
  fives_.y() = 5;
  delta_.x() = 10;
  delta_.y() = 10;
  DEMfile_ = "USGS_A15_Q111_LRO_NAC_DEM_26N004E_150cmp.30mpp.tif";
  }

  virtual void TearDown() {}

  string DEMfile_;
  Vector2 fives_;
  Vector2 delta_;
};


TEST_F( Features_Test, GetFeatures ){
  boost::shared_ptr<DiskImageResource> fore_rsrc( new DiskImageResourceGDAL( DEMfile_ ) );
  float noDataVal = fore_rsrc->nodata_read();
  DiskImageView<PixelGray<float> > foreDEM( fore_rsrc );
  GeoReference foreDEMGeo;
  read_georeference( foreDEMGeo, DEMfile_ );

  vector<Vector3> features = GetFeatures( foreDEM, foreDEMGeo, foreDEM, foreDEMGeo, 
                                          fives_, delta_, noDataVal );

  ASSERT_EQ( (unsigned int)5428, features.size()  ) << "There is an incorrect number of features.";
  EXPECT_NEAR( 1.35555e+06, features[100].x(), 1) << "Test feature has wrong x.";
  EXPECT_NEAR( 328576,      features[100].y(), 1) << "Test feature has wrong y.";
  EXPECT_NEAR( 1.0327e+06,  features[100].z(), 5) << "Test feature has wrong z.";
}

TEST_F( Features_Test, FindMatches ){
  boost::shared_ptr<DiskImageResource> fore_rsrc( new DiskImageResourceGDAL( DEMfile_ ) );
  float noDataVal = fore_rsrc->nodata_read();
  DiskImageView<PixelGray<float> > foreDEM( fore_rsrc );
  GeoReference foreDEMGeo;
  read_georeference( foreDEMGeo, DEMfile_ );


  Vector2 delta( 0, 0 );
  vector<Vector3> features = GetFeatures( foreDEM, foreDEMGeo,
                                          fives_, delta, noDataVal );

  vector<Vector3> matches( features.size() );
 
  Vector2 window( 5, 5 );
  FindMatches( features, foreDEM, foreDEMGeo, foreDEMGeo, matches, window, noDataVal );

  int diff_count = 0;
  for( unsigned int i=0; i < matches.size(); ++i ){
    if( !(matches[i] == features[i]) ){
      cerr << matches[i] << " " << features[i] << endl;
      ++diff_count;
    }
  }
  ASSERT_EQ( 0, diff_count ) << "Matched features differ from reference.";
}

TEST( ICP_DEM_2_DEM, works ) {
  string DEMfile( "USGS_A15_Q111_LRO_NAC_DEM_26N004E_150cmp.30mpp.tif" );
  boost::shared_ptr<DiskImageResource> fore_rsrc( new DiskImageResourceGDAL( DEMfile ) );
  float noDataVal = fore_rsrc->nodata_read();
  DiskImageView<PixelGray<float> > foreDEM( fore_rsrc );
  GeoReference foreDEMGeo;
  read_georeference( foreDEMGeo, DEMfile );

  vw_log().console_log().rule_set().add_rule(vw::InfoMessage,"*");

  Vector2 delta( 0, 0 );
  Vector2 fives( 5, 5 );
  vector<Vector3> features = GetFeatures( foreDEM, foreDEMGeo,
                                          fives, delta, noDataVal );
  for( unsigned int i = 0; i < features.size(); ++i){
     features[i].x() += 5;
  }
  
  Vector3 translation;
  Matrix<float, 3,3 > rotation = identity_matrix(3);
  Vector3 center = find_centroid( features );
  //vector<float> errors( features.size() );
  float error;

  struct CoregistrationParams settings;
  settings.maxNumIter = 10;
  settings.minConvThresh = 0.01;
  settings.matchWindowHalfSize = fives;
  settings.noDataVal = noDataVal;

  ICP_DEM_2_DEM( features, foreDEM, foreDEMGeo, foreDEMGeo, 3,
                 settings, translation, rotation, center, error );
  cout << translation << endl;
  EXPECT_NEAR( -5.0, translation.x(), 0.1 ) << "Translation in X is wrong.";
  EXPECT_NEAR( 0.0, translation.y(), 0.1 ) << "Translation in Y is wrong.";
  EXPECT_NEAR( 0.0, translation.z(), 0.1 ) << "Translation in Z is wrong.";
  for( unsigned int i = 0; i < 3; ++i ){
    for( unsigned int j = 0; j < 3; j++ ){
      float truth = 0.0;
      if( i == j ){ truth = 1.0; }
      EXPECT_NEAR( truth, rotation(i,j), 0.001 ) << "Rotation matrix is wrong.";
    }
  }
  EXPECT_LE( error, settings.minConvThresh ) << "The error value is wrong.";
}

class DEM_Matching : public ::testing::Test {
  protected:
  virtual void SetUp(){
    DEMfile_ = "USGS_A15_Q111_LRO_NAC_DEM_26N004E_150cmp.30mpp.tif";
    read_georeference( DEMGeo_, DEMfile_ );
    shots_ = LOLAFileRead( "RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv" );
  }

  virtual void TearDown() {}

  string DEMfile_;
  GeoReference DEMGeo_;
  vector<vector<LOLAShot> > shots_;
};

TEST_F( DEM_Matching, FindMatchesFromDEM ){
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL( DEMfile_ ) );
  double nodata = rsrc->nodata_read();
  DiskImageView<float> DEM(rsrc);

  GetAllPtsFromDEM( shots_, DEM, DEMGeo_, nodata );

  vector<Vector3> xyz;//LOLA
  vector<Vector3> llr;
  for( unsigned int i = 0; i < shots_.size(); ++i ){
    for( unsigned int j = 0; j < shots_[i].size(); ++j ){
    if( shots_[i][j].DEMPt[0].valid ){
      Vector3 model;
      model = shots_[i][j].LOLAPt[2];
	  model.z() *= 1000; // LOLA data is in km, DEMGeo is in m (for LROC DTMs).
      Vector3 xyzModel = lon_lat_radius_to_xyz(model);
      xyz.push_back(xyzModel);
      llr.push_back(model);
      }
    }
  }

  vector<Vector3> matches( xyz.size() );
  Vector3 modelCentroid = find_centroid( xyz );
  Matrix<float, 3,3 > rotation = identity_matrix(3);
  Vector3 translation;
  Vector2 window( 5, 5 );

  FindMatchesFromDEM( xyz, llr, DEM, DEMGeo_, matches, translation, 
                      rotation, modelCentroid, nodata, window );

  ASSERT_EQ( xyz.size(), matches.size() ) << "The match vector is the wrong size.";
  EXPECT_NEAR( matches[761].x(), 1.56342e+06, 10 ) << "The x value is off.";
  EXPECT_NEAR( matches[761].y(), 98450.7, .1 ) << "The y value is off.";
  EXPECT_NEAR( matches[761].z(), 752184, 1 ) << "The y value is off.";
  //for( unsigned int i = 0; i < xyz.size(); ++i ){
  //  if( xyz[i] != matches[i] ) { cerr << i << " " << xyz[i] << " " << matches[i] << endl; }
  //}
}

TEST_F( DEM_Matching, ICP_LIDAR_2_DEM ){
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL( DEMfile_ ) );
  double nodata = rsrc->nodata_read();
  DiskImageView<float> DEM(rsrc);

  GetAllPtsFromDEM( shots_, DEM, DEMGeo_, nodata );

  vector<Vector3> xyz;//LOLA
  vector<Vector3> llr;
  for( unsigned int i = 0; i < shots_.size(); ++i ){
    for( unsigned int j = 0; j < shots_[i].size(); ++j ){
    if( shots_[i][j].DEMPt[0].valid ){
      Vector3 model;
      model = shots_[i][j].LOLAPt[2];
	  model.z() *= 1000; // LOLA data is in km, DEMGeo is in m (for LROC DTMs).
      Vector3 xyzModel = lon_lat_radius_to_xyz(model);
      xyz.push_back(xyzModel);
      llr.push_back(model);
      }
    }
  }

  vector<Vector3> xyzMatchArray( xyz.size() );
  valarray<float> xyzErrorArray( xyz.size() );
  Matrix<float, 3,3 > rotation = identity_matrix(3);
  Vector3 translation;
  Vector3 modelCentroid = find_centroid( xyz );
  struct CoregistrationParams settings;
  settings.maxNumIter = 10;
  settings.minConvThresh = 0.01;
  settings.matchWindowHalfSize = Vector2( 5, 5 );
  settings.noDataVal = nodata;

  //vw_log().console_log().rule_set().add_rule(vw::InfoMessage,"*");

  ICP_LIDAR_2_DEM( xyzMatchArray, DEM, DEMGeo_, xyz, llr, settings, 
                   translation, rotation, modelCentroid, xyzErrorArray);

  EXPECT_NEAR( translation.x(), -2.87517, 0.01 ) << "Translation in X is wrong.";
  EXPECT_NEAR( translation.y(), -10.2501, 0.01 ) << "Translation in Y is wrong.";
  EXPECT_NEAR( translation.z(), -0.36370, 0.01 ) << "Translation in Z is wrong.";
  for( unsigned int i = 0; i < 3; ++i ){
    for( unsigned int j = 0; j < 3; j++ ){
      float truth = 0.0;
      if( i == j ){ truth = 1.0; }
      EXPECT_NEAR( truth, rotation(i,j), 0.01 ) << "Rotation matrix is wrong.";
    }
  }
  EXPECT_LE( xyzErrorArray.sum()/xyzErrorArray.size(), 14.5 ) << "The error value is wrong.";
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
