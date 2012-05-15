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

TEST( ComputeDEMRotation, works ){
  vector<Vector3> features;
  features.push_back( Vector3( 5,5,6 ) );
  features.push_back( Vector3( 5,6,5 ) );
  features.push_back( Vector3( 6,5,5 ) );
  vector<Vector3> match;
  match.push_back( Vector3( 5,5,4 ) );
  match.push_back( Vector3( 5,6,5 ) );
  match.push_back( Vector3( 6,5,5 ) );

  Matrix<float, 3, 3> truth;
  truth(0,0) =  1.0/3.0;
  truth(0,1) = -2.0/3.0;
  truth(0,2) = -2.0/3.0;
  truth(1,0) = -2.0/3.0;
  truth(1,1) =  1.0/3.0;
  truth(1,2) = -2.0/3.0;
  truth(2,0) =  2.0/3.0;
  truth(2,1) =  2.0/3.0;
  truth(2,2) = -1.0/3.0;
  
  Vector3 f_center = find_centroid( features );
  Vector3 m_center = find_centroid( match );
  Matrix<float, 3, 3> test = ComputeDEMRotation( features, match, f_center, m_center );

  for( unsigned int i = 0; i < 3; ++i ){
    for( unsigned int j = 0; j < 3; ++j ){
    EXPECT_NEAR( truth(i,j), test(i,j), 0.0001 ) << "Rotation elements aren't correct.";
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
