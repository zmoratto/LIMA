#include <boost/filesystem.hpp>

#include "../tracks.h"
#include "gtest/gtest.h"

// #include <asp/IsisIO.h>
// #include <asp/IsisIO/IsisCameraModel.h>

namespace fs = boost::filesystem;

using namespace std;

TEST( LOLAfileRead_Test, canRead ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<unsigned int> truth_sizes(5);
  truth_sizes[0] = 422;
  truth_sizes[1] = 1299;
  truth_sizes[2] = 1431;
  truth_sizes[3] = 1568;
  truth_sizes[4] = 1325;

  vector<vector<LOLAShot> > test = CSVFileRead( p.string() );

  ASSERT_EQ( (unsigned int)5, test.size() ) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth_sizes.size(); ++i ) {
    EXPECT_EQ(truth_sizes[i], test[i].size()) << "Vectors truth and test differ at index " << i;
  }

  // Select one point to test against.
  ASSERT_EQ( (unsigned int)5, test[2][0].LOLAPt.size() )<< "Shot doesn't have enough points.";

  ASSERT_EQ( 3, test[2][0].LOLAPt[2].s ) << "This point doesn't have the right detector.";
  ASSERT_EQ( 2009, test[2][0].LOLAPt[2].year ) << "This point doesn't have the right year.";
  ASSERT_EQ( 10, test[2][0].LOLAPt[2].month ) << "This point doesn't have the right month.";
  ASSERT_EQ( 3, test[2][0].LOLAPt[2].day) << "This point doesn't have the right day.";
  ASSERT_EQ( 13, test[2][0].LOLAPt[2].hour) << "This point doesn't have the right hour.";
  ASSERT_EQ( 26, test[2][0].LOLAPt[2].min) << "This point doesn't have the right minute.";
  ASSERT_NEAR( 45.2282, test[2][0].LOLAPt[2].sec, 0.0001) << "This point doesn't have the right seconds.";
  
}

TEST( GetPointFromIndex_Test, works ){
  pointCloud point0( Vector3(1,1,1), 2009, 10, 3, 13, 26, 0, 0 );
  pointCloud point1( Vector3(1,1,1), 2009, 10, 3, 13, 26, 0, 1 );
  pointCloud point2( Vector3(1,1,1), 2009, 10, 3, 13, 26, 0, 2 );
  pointCloud point3( Vector3(1,1,1), 2009, 10, 3, 13, 26, 0, 3 );
  pointCloud point4( Vector3(1,1,1), 2009, 10, 3, 13, 26, 0, 4 );
  vector<pointCloud> truth(5);
  truth[0] = point0;
  truth[1] = point1;
  truth[2] = point2;
  truth[3] = point3;
  truth[4] = point4;

  pointCloud test = GetPointFromIndex( truth, 3 );

  ASSERT_EQ( point3.s, test.s ) << "Point with the wrong detector.";
}

TEST( FindMinMaxLat_Test, works ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<vector<LOLAShot> > test = LOLAFileRead( p.string() );

  Vector4 bounds = FindMinMaxLat( test );
  //Vector4(24.5005,27.5,3.06477,3.77167)
  EXPECT_NEAR( 24.5005, bounds[0], 0.0001 ) << "Minimum latitude is wrong.";
  EXPECT_NEAR( 27.5   , bounds[1], 0.0001 ) << "Maximum latitude is wrong.";
  EXPECT_NEAR( 3.06477, bounds[2], 0.0001 ) << "Minimum longitude is wrong.";
  EXPECT_NEAR( 3.77167, bounds[3], 0.0001 ) << "Maximum longitude is wrong.";
}

TEST( ComputeAllReflectance_Test, works ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<vector<LOLAShot> > shots = LOLAFileRead( p.string() );
  unsigned int num_shots = shots.size();

  /*
  camera::IsisCameraModel model("M111578606RE.10mpp.cub");
  Vector2 pixel_location( model.samples()/2, model.lines()/2 );
  Vector3 cameraPosition = model.camera_center( pixel_location );
  Vector3 lightPosition = model.sun_position( pixel_location );
  
  cout << "camera " << cameraPosition << endl;
  // camera Vector3(1.60109e+06,82708.6,784432)
  cout << "light " << lightPosition << endl;
  // light Vector3(1.28119e+11,7.5658e+10,-3.98442e+09)
  */

  Vector3 cameraPosition( 1.60109e+06,82708.6,784432 );
  Vector3 lightPosition( 1.28119e+11,7.5658e+10,-3.98442e+09 );
  struct CoregistrationParams settings;

  int test_valid_points = ComputeAllReflectance( shots,  
                                                 cameraPosition, 
                                                 lightPosition, 
                                                 settings);
  ASSERT_EQ( num_shots, shots.size() ) << "The length of the vector of shots was altered.";
  ASSERT_EQ( 2723, test_valid_points ) << "There is the wrong number of valid points.";
  ASSERT_NEAR( 0.859062, shots[4][1306].reflectance, 0.00001 ) << "The reflectance is wrong.";
}

TEST( ComputeGainBiasFactor, DISABLED_vector_of_shots ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<vector<LOLAShot> > shots = LOLAFileRead( p.string() );
  // I think we need to run ComputeAllReflectance here, but it needs a test, first.
  Vector2 test = ComputeGainBiasFactor( shots );

  cout << test << endl;
}

TEST( ComputeGainBiasFactor, DISABLED_vector_of_vector_of_shots ){
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
