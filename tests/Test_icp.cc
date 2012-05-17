#include <boost/filesystem.hpp>

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


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
