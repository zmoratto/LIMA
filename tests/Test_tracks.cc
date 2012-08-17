#include <boost/filesystem.hpp>

#include "../tracks.h"
#include "../match.h" // For FindMatches2D()
#include "../featuresLOLA.h" // For ComputeSalientLOLAFeature()
#include "gtest/gtest.h"

// #include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

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

TEST( GetAllPtsFromImage, works ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<vector<LOLAShot> > shots = LOLAFileRead( p.string() );
  unsigned int num_shots = shots.size();
  fs::path im("M111578606RE.10mpp.tif");

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL( im.string() ) );
  DiskImageView<PixelGray<uint8> > DRG( rsrc );
  GeoReference DRGGeo;
  read_georeference(DRGGeo, im.string() );

  GetAllPtsFromImage( shots, DRG, DRGGeo );
  ASSERT_EQ( num_shots, shots.size() ) << "The length of the vector of shots was altered.";

  int valid_counter = 0;
  for( unsigned int i = 0; i < shots.size(); ++i ){
    for( unsigned int j = 0; j < shots[i].size(); ++j ){
      if( shots[i][j].valid == 1 ){ 
        ++valid_counter;
        // for( unsigned int k = 0; k < shots[i][j].imgPt.size(); ++k ){
        //   cout << i << " " << j << " " << k
        //        << " x: " << shots[i][j].imgPt[k].x 
        //        << " y: " << shots[i][j].imgPt[k].y
        //        << " val: " << shots[i][j].imgPt[k].val << endl;
      }
    }
  }
  ASSERT_EQ( 331, valid_counter ) << "There is a different number of valid shots.";
  ASSERT_EQ( 1, shots[4][785].valid ) << "Shot isn't valid.";
  ASSERT_EQ( (unsigned int)5, shots[4][785].imgPt.size() ) << "Vector of imgPt wrong size.";
  ASSERT_NEAR( 234.945, shots[4][785].imgPt[3].x, 0.001 ) << "The x value is wrong.";
  ASSERT_NEAR( 1826.62, shots[4][785].imgPt[3].y, 0.01 ) << "The y value is wrong.";
  ASSERT_NEAR( 10, shots[4][785].imgPt[3].val, 0.1 ) << "The y value is wrong.";

  //SaveImagePoints( shots, 1, "SaveImagePoints_Test.txt" );
  //SaveAltitudePoints( shots, 1, "SaveAltitudePoints_Test.txt" );
  //SaveDEMPoints( shots, "M111578606RE.10mpp.tif" , "SaveDEMPoints_Test.txt" );
}

TEST( GetAllPtsFromCub, works ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<vector<LOLAShot> > shots = LOLAFileRead( p.string() );
  unsigned int num_shots = shots.size();
  fs::path cub("M111578606RE.10mpp.cub");

  boost::shared_ptr<DiskImageResource> rsrc(new DiskImageResourceIsis(cub.string()));
  DiskImageView<PixelGray<float> > cubView(rsrc);
  ImageView<PixelGray<float> > cubImage(cubView.cols(), cubView.rows());
  cubImage = cubView;
  camera::IsisCameraModel model(cub.string());

  GetAllPtsFromCub( shots, model, cubImage);
  ASSERT_EQ( num_shots, shots.size() ) << "The length of the vector of shots was altered.";

  int valid_counter = 0;
  for( unsigned int i = 0; i < shots.size(); ++i ){
    for( unsigned int j = 0; j < shots[i].size(); ++j ){
      if( shots[i][j].valid == 1 ){ 
        ++valid_counter;
        // for( unsigned int k = 0; k < shots[i][j].imgPt.size(); ++k ){
        //   cout << i << " " << j << " " << k
        //        << " x: " << shots[i][j].imgPt[k].x 
        //        << " y: " << shots[i][j].imgPt[k].y
        //        << " val: " << shots[i][j].imgPt[k].val << endl;
      }
    }
  }
  ASSERT_EQ( 301, valid_counter ) << "There is a different number of valid shots.";
  ASSERT_EQ( 1, shots[4][785].valid ) << "Shot isn't valid.";
  ASSERT_EQ( (unsigned int)5, shots[4][785].imgPt.size() ) << "Vector of imgPt wrong size.";
  ASSERT_NEAR( 222.973, shots[4][785].imgPt[3].x, 0.001 ) << "The x value is wrong.";
  ASSERT_NEAR( 1826.39, shots[4][785].imgPt[3].y, 0.01 ) << "The y value is wrong.";
  ASSERT_NEAR( 0.0439, shots[4][785].imgPt[3].val, 0.0001 ) << "The image value is wrong.";
}


class ReflectanceTests : public ::testing::Test {
  protected:
  virtual void SetUp() {
    p_ = "RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv";
    camera_pos_(0) = 1.60109e+06;
    camera_pos_(1) = 82708.6;
    camera_pos_(2) = 784432;
    light_pos_(0) = 1.28119e+11;
    light_pos_(1) = 7.5658e+10;
    light_pos_(2) = -3.98442e+09;
  }

  virtual void TearDown() {}

  fs::path p_;
  Vector3 camera_pos_;
  Vector3 light_pos_;
};

TEST_F( ReflectanceTests, ComputeAllReflectance_Test ){
  vector<vector<LOLAShot> > shots = LOLAFileRead( p_.string() );
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

  int test_valid_points = ComputeAllReflectance( shots,  
                                                 camera_pos_, 
                                                 light_pos_);
  ASSERT_EQ( num_shots, shots.size() ) << "The length of the vector of shots was altered.";
  ASSERT_EQ( 2787, test_valid_points ) << "There is the wrong number of valid points.";
  ASSERT_NEAR( 0.87029075622558594, shots[4][1306].reflectance, 0.00001 ) << "The reflectance is wrong.";
}

TEST_F( ReflectanceTests,  ComputeGainBiasFactor_vector_of_shots ){
  vector<vector<LOLAShot> > shots = LOLAFileRead( p_.string() );

  fs::path im("M111578606RE.10mpp.tif");
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL( im.string() ) );
  DiskImageView<PixelGray<uint8> > DRG( rsrc );
  GeoReference DRGGeo;
  read_georeference(DRGGeo, im.string() );
  GetAllPtsFromImage( shots, DRG, DRGGeo );

  ComputeAllReflectance( shots, camera_pos_, light_pos_);
  vector<vector< AlignedLOLAShot> > aligned = initialize_aligned_lola_shots(shots);
  Vector2 test = ComputeGainBiasFactor( aligned[4] );

  ASSERT_NEAR( 11.069457054138184, test[0], 0.00001 ) << "Gain is different.";
  ASSERT_NEAR( 0.0, test[1], 0.0001 ) << "Bias is different.";
}

TEST_F( ReflectanceTests, ComputeGainBiasFactor_vector_of_vector_of_shots ){
  vector<vector<LOLAShot> > shots = LOLAFileRead( p_.string() );
  ComputeAllReflectance( shots,  camera_pos_, light_pos_ );

  fs::path im("M111578606RE.10mpp.tif");
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL( im.string() ) );
  DiskImageView<PixelGray<uint8> > DRG( rsrc );
  GeoReference DRGGeo;
  read_georeference(DRGGeo, im.string() );
  GetAllPtsFromImage( shots, DRG, DRGGeo );

  // Need to add some other reflectances.  Only one track seems to have valid
  // values coming out of ComputeAllReflectance(), so we will fake some.
  for( unsigned int j = 0; j < shots[2].size(); ++j ){
    shots[2][j].valid = 1;
    shots[2][j].reflectance = 0.3;
  }

  vector<vector< AlignedLOLAShot> > aligned = initialize_aligned_lola_shots(shots);
  Vector2 test = ComputeGainBiasFactor( aligned );

  EXPECT_NEAR( 4.2450780868530273,  test[0], 0.0001 ) << "Gain is different.";
  EXPECT_NEAR( 0.0, test[1], 0.0001 ) << "Bias is different.";

  //SaveReflectancePoints( shots, test, "SaveReflectancePoints_Test.txt");
}

TEST( ComputeMinMaxValuesFromCub, works ){
  boost::shared_ptr<DiskImageResource> rsrc(new DiskImageResourceIsis("AS15-M-2327.lev1.500.cub"));
  DiskImageView<PixelGray<float> > cubView(rsrc);
  ImageView<PixelGray<float> > cubImage(cubView.cols(), cubView.rows());
  cubImage = cubView;
  Vector2 test = ComputeMinMaxValuesFromCub(cubImage );

  // The output from the ISIS stats program is:
  // Minimum                 = -1.29148247651756e-09
  // Maximum                 = 65535.000000001
  //
  // The output from this function is
  // min=-32752, max=32767
  // so this function is not getting the "ISIS values" but the raw values.

  EXPECT_EQ( -32752, test[0] ) << "Wrong min.";
  EXPECT_EQ(  32767, test[1] ) << "Wrong max.";

}

class GetAllPtsFromDEM_Test : public ::testing::Test {
  protected:
  virtual void SetUp() {
    p_ = "RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv";
    d_ = "USGS_A15_Q111_LRO_NAC_DEM_26N004E_150cmp.30mpp.tif";
    nodata_ = -10000;
  }

  virtual void TearDown() {}

  fs::path p_;
  fs::path d_;
  double nodata_;
};


TEST_F( GetAllPtsFromDEM_Test, works ){
  vector<vector<LOLAShot> > shots = LOLAFileRead( p_.string() );

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(d_.string()) );
  if( rsrc->has_nodata_read() ){ nodata_ = rsrc->nodata_read(); }
  DiskImageView<float> DEM(rsrc);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, d_.string() );

  GetAllPtsFromDEM( shots, DEM, DEMGeo, nodata_ );
  
  int valid_counter = 0;
  for( unsigned int i = 0; i < shots.size(); ++i ){
    for( unsigned int j = 0; j < shots[i].size(); ++j ){
      if( shots[i][j].valid == 1 ){ 
        for( unsigned int k = 0; k < shots[i][j].DEMPt.size(); ++k ){
          if( shots[i][j].DEMPt[k].valid ){
          ++valid_counter;
          //cout << i << " " << j << " " << k
          //     << " x: " << shots[i][j].DEMPt[k].x 
          //     << " y: " << shots[i][j].DEMPt[k].y
          //     << " val: " << shots[i][j].DEMPt[k].val 
          //     << " valid: " << shots[i][j].DEMPt[k].valid << endl;
          }
        }
      }
    }
  }

  ASSERT_EQ( 3927, valid_counter ) << "There is a different number of valid DEMPts.";
  ASSERT_EQ( 1, shots[2][762].DEMPt[3].valid ) << "Shot isn't valid.";
  ASSERT_EQ( (unsigned int)5, shots[2][762].DEMPt.size() ) << "Vector of imgPt wrong size.";
  ASSERT_NEAR( 110.213, shots[2][762].DEMPt[3].x, 0.001 ) << "The x value is wrong.";
  ASSERT_NEAR( 774.582, shots[2][762].DEMPt[3].y, 0.001 ) << "The y value is wrong.";
  ASSERT_NEAR( 1738.1,  shots[2][762].DEMPt[3].val, 0.1 ) << "The elevation value is wrong.";
}

TEST_F( GetAllPtsFromDEM_Test, precision ){
  vector<vector<LOLAShot> > shots = LOLAFileRead( p_.string() );

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(d_.string()) );
  if( rsrc->has_nodata_read() ){ nodata_ = rsrc->nodata_read(); }
  DiskImageView<float> DEM(rsrc);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, d_.string() );

  boost::shared_ptr<DiskImageResource> p_rsrc( new DiskImageResourceGDAL("USGS_A15_Q111_LRO_NAC_DEM_26N004E_150cmp.30mpp.pslope.tif") );
  DiskImageView<float> DEM_Prec(p_rsrc);
  

  GetAllPtsFromDEM_Prec( shots, DEM, DEMGeo, nodata_, DEM_Prec );
  
  int valid_counter = 0;
  for( unsigned int i = 0; i < shots.size(); ++i ){
    for( unsigned int j = 0; j < shots[i].size(); ++j ){
      if( shots[i][j].valid == 1 ){ 
        for( unsigned int k = 0; k < shots[i][j].DEMPt.size(); ++k ){
          if( shots[i][j].DEMPt[k].valid ){
          ++valid_counter;
          // cout << i << " " << j << " " << k
          //      << " x: " << shots[i][j].DEMPt[k].x 
          //      << " y: " << shots[i][j].DEMPt[k].y
          //      << " val: " << shots[i][j].DEMPt[k].val 
          //      << " valid: " << shots[i][j].DEMPt[k].valid << endl;
          }
        }
      }
    }
  }

  ASSERT_EQ( 2473, valid_counter ) << "There is a different number of valid DEMPts.";
  ASSERT_EQ( 1, shots[4][523].DEMPt[3].valid ) << "Shot isn't valid.";
  ASSERT_EQ( (unsigned int)5, shots[4][523].DEMPt.size() ) << "Vector of imgPt wrong size.";
  ASSERT_NEAR( 80.3372, shots[4][523].DEMPt[3].x, 0.001 ) << "The x value is wrong.";
  ASSERT_NEAR( 32.4886, shots[4][523].DEMPt[3].y, 0.001 ) << "The y value is wrong.";
  ASSERT_NEAR( 1735.5,  shots[4][523].DEMPt[3].val, 0.1 ) << "The elevation value is wrong.";
}

TEST( GCP, Update ){
  fs::path p("RDR_3E4E_24N27NPointPerRow_csv_table-truncated.csv");
  vector<vector<LOLAShot> > shots = LOLAFileRead( p.string() );

  vector<float> filter = MakeLidarFilter( 12 );

  for( unsigned int i = 0; i < shots.size(); ++i ){
    ComputeSalientLOLAFeature( shots[i], filter, 0.008 );
  }

  vector<gcp> gcpArray;
  for (unsigned int t=0; t < shots.size(); t++){
    for (unsigned int s=0; s < shots[t].size(); s++){
      if (shots[t][s].featurePtLOLA==1){
        gcp this_gcp;
        this_gcp.lon = shots[t][s].LOLAPt[2].x();
        this_gcp.lat = shots[t][s].LOLAPt[2].y(); 
        this_gcp.rad = shots[t][s].LOLAPt[2].z()*1000;
        this_gcp.sigma_lon = 1.0;
        this_gcp.sigma_lat = 1.0;
        this_gcp.sigma_rad = 1.0;
        gcpArray.push_back(this_gcp);
      }
    }
  }

  string map( "M111578606RE.10mpp.cub" );
  string cam( "cam_cub_filename" );
  
  boost::shared_ptr<DiskImageResource> rsrc(new DiskImageResourceIsis(map));
  DiskImageView<PixelGray<float> > cubView(rsrc);
  ImageView<PixelGray<float> > cubImage(cubView.cols(), cubView.rows());
  cubImage = cubView;
  camera::IsisCameraModel model(map);


  
  GetAllPtsFromCub( shots, model, cubImage);

  Vector2 origin;
  Vector3 cameraPosition = model.camera_center( origin );
  Vector3 lightPosition = model.sun_position( origin );
  ComputeAllReflectance( shots,  cameraPosition, lightPosition); // use bogus positions?
  
  vector<float> errors;
  vector<vector< AlignedLOLAShot> > aligned = initialize_aligned_lola_shots(shots);
  vector<Vector4> matches = FindMatches2D( aligned, map, 5, 80, errors );

  UpdateGCP( shots, matches, cam, map, gcpArray, 4.0 ); 

  ASSERT_EQ( (unsigned int)23, gcpArray.size() ) << "Array of GCPs not the right size.";
  EXPECT_EQ( cam, gcpArray[13].filename[0] ) << "Wrong camera file name.";
  EXPECT_NEAR( 1609.619, gcpArray[14].x[0], 0.001 ) << "Wrong x value.";
  EXPECT_NEAR( 8.0933, gcpArray[14].y[0], 0.01 ) << "Wrong y value.";
  EXPECT_NEAR( 527.185, gcpArray[14].x_before[0], 0.001 ) << "Wrong x_before value.";
  EXPECT_NEAR( 7075.3, gcpArray[14].y_before[0], 0.1 ) << "Wrong y_before value.";
  EXPECT_EQ( 4, gcpArray[14].trackIndex ) << "Wrong trackIndex value.";
  EXPECT_EQ( 747, gcpArray[14].shotIndex ) << "Wrong shotIndex value.";

  // string saveGCP( "SaveGCPoints_test.txt" );
  // SaveGCPoints( gcpArray, saveGCP );
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
