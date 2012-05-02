#include <boost/filesystem.hpp>

#include "../util.h"
#include "gtest/gtest.h"

namespace fs = boost::filesystem;

using namespace std;

void test_write_vector_to_stream( const vector<string>& v, std::ostream& stream ) {
  for( vector<string>::const_iterator i = v.begin(); 
       i != v.end();
       ++i ){
    stream << *i << endl;
  }
  return;
}

class ReadFileToVectorTest : public ::testing::Test {
  protected:
  virtual void SetUp() {
    truth_.resize(4);
    truth_[0] = "This";
    truth_[1] = "is";
    truth_[2] = "a";
    truth_[3] = "test.";

    truthi_.resize(4);
    truthi_[0] = 1;
    truthi_[1] = 1;
    truthi_[2] = 2;
    truthi_[3] = 3;
    

    p_ = "ReadFileToVector_temp.txt";
    pi_ = "ReadFileToVector_tempi.txt";

    std::ofstream file( p_.c_str() );
    ASSERT_TRUE( file ) << "Can't open file " << p_ ;
    test_write_vector_to_stream( truth_, file );
    file.close();

    std::ofstream filei( pi_.c_str() );
    ASSERT_TRUE( filei ) << "Can't open file " << pi_ ;
    for( vector<int>::const_iterator i = truthi_.begin(); i != truthi_.end(); ++i ) {
      filei << *i << endl;
    }
    filei.close();
  }

  virtual void TearDown() {
    fs::remove( p_ );
    fs::remove( pi_ );
  }

  vector<string> truth_;
  vector<int> truthi_;
  fs::path p_;
  fs::path pi_;
};

TEST_F( ReadFileToVectorTest, stringvector_returnsString ) {
  vector<string> test = ReadVectorFrom<string>( p_.string() );

  ASSERT_EQ(truth_.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth_.size(); ++i ) {
    EXPECT_EQ(truth_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( ReadFileToVectorTest, stringvector_returnsVoid ) {
  vector<string> test;
  ReadFileList( p_.string(), test );

  ASSERT_EQ(truth_.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth_.size(); ++i ) {
    EXPECT_EQ(truth_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( ReadFileToVectorTest, intvector_returnsint ) {
  vector<int> test;
  int return_value = ReadOverlapList( pi_.string(), test );

  ASSERT_EQ( 1, return_value ) << "Function could not open file " << pi_.string();

  ASSERT_EQ(truthi_.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth_.size(); ++i ) {
    EXPECT_EQ(truthi_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}



class AccessDataFilesFromInputTest: public ::testing::Test {
  protected:
  virtual void SetUp() {
    truth1_.push_back( "File1" );

    truth2_.resize(3);
    truth2_[0] = "File01"; 
    truth2_[1] = "File02"; 
    truth2_[2] = "File03.data";

    truth1_p_ = "AccessDataFilesFromInput_truth1.txt";
    truth2_p_ = "AccessDataFilesFromInput_truth2.txt";

    list1_.push_back( truth1_p_.string() );

    list2_.resize(4);
    list2_[0] = "File_A"; 
    list2_[1] = "File_B"; 
    list2_[2] = "File_C.data"; 
    list2_[3] = truth2_p_.string();

    truth3_ = list2_;
    truth3_.pop_back();
    truth3_.insert( truth3_.end(), truth2_.begin(), truth2_.end() );

    // Write out some files to test the file reading part
    std::ofstream file0( truth1_p_.c_str() );
    ASSERT_TRUE( file0 ) << "Can't open file " << truth1_p_ ;
    test_write_vector_to_stream( truth1_, file0 );
    file0.close();

    std::ofstream file1( truth2_p_.c_str() );
    ASSERT_TRUE( file1 ) << "Can't open file " << truth2_p_ ;
    test_write_vector_to_stream( truth2_, file1 );
    file1.close();
  }

  virtual void TearDown() {
    fs::remove( truth1_p_ );
    fs::remove( truth2_p_ );
  }

  vector<string> truth1_;
  vector<string> truth2_;
  vector<string> truth3_;

  vector<string> list1_;
  vector<string> list2_;

  fs::path truth1_p_;
  fs::path truth2_p_;
  
};

TEST_F( AccessDataFilesFromInputTest, singleElement ) {
  vector<string> test;
  test = AccessDataFilesFromInput( truth1_ );
  ASSERT_EQ(truth1_.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth1_.size(); ++i ) {
    EXPECT_EQ(truth1_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( AccessDataFilesFromInputTest, singleElementAsString ) {
  vector<string> test;
  test = AccessDataFilesFromInput( truth1_.front() );
  ASSERT_EQ(truth1_.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth1_.size(); ++i ) {
    EXPECT_EQ(truth1_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( AccessDataFilesFromInputTest, multipleElements ) {
  vector<string> test;
  test = AccessDataFilesFromInput( truth2_ );
  ASSERT_EQ(truth2_.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth2_.size(); ++i ) {
    EXPECT_EQ(truth2_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( AccessDataFilesFromInputTest, singleFileName ) {
  vector<string> test;
  test = AccessDataFilesFromInput( list1_ );
  ASSERT_EQ(truth1_.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth1_.size(); ++i ) {
    EXPECT_EQ(truth1_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( AccessDataFilesFromInputTest, singleFileNameAsString ) {
  vector<string> test;
  test = AccessDataFilesFromInput( list1_.front() );
  ASSERT_EQ(truth1_.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth1_.size(); ++i ) {
    EXPECT_EQ(truth1_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( AccessDataFilesFromInputTest, multipleElementsWFileName ) {
  vector<string> test;
  test = AccessDataFilesFromInput( list2_ );
  ASSERT_EQ(truth3_.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth3_.size(); ++i ) {
    EXPECT_EQ(truth3_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST( ComputeGeoBoundary_Test, temporary ) {
  string filename( "AS15-M-2327.lev1.500.cub" );
  Vector4 test = ComputeGeoBoundary( filename );
  unsigned int four = 4;
  ASSERT_EQ(four, test.size()) << "Vectors are of unequal length";

  float minLon = test(0);
  float maxLon = test(1);
  float minLat = test(2);
  float maxLat = test(3);

  //  cout << test << endl;
  // Vector4(-31.2156,-24.3621,25.4635,31.4475)
  EXPECT_NEAR(-31.2156,minLon,0.0001) << "Minimum Longitude is incorrect."; 
  EXPECT_NEAR(-24.3621,maxLon,0.0001) << "Maximum Longitude is incorrect."; 
  EXPECT_NEAR( 25.4635,minLat,0.0001) << "Minimum Latitude is incorrect."; 
  EXPECT_NEAR( 31.4475,maxLat,0.0001) << "Maximum Latitude is incorrect."; 
}

class makeOverlapList_Test : public ::testing::Test {
  protected:
  virtual void SetUp() {
    files_.push_back( "AS15-M-2327.lev1.500.cub" );
    files_.push_back( "AS15-M-2327.lev1.500.cub" );
    // Vector4(-31.2156,-24.3621,25.4635,31.4475)
    inside_(0) = -30;
    inside_(1) = -25;
    inside_(2) =  26;
    inside_(3) =  30;
    partial_(0) = -28;
    partial_(1) = -20;
    partial_(2) =  27;
    partial_(3) =  29;
    outside_(0) = 0;
    outside_(1) = 1;
    outside_(2) = 0;
    outside_(3) = 1;
  }

  virtual void TearDown() {}

  vector<string> files_;
  Vector4 inside_;
  Vector4 partial_;
  Vector4 outside_;
};
    

TEST_F( makeOverlapList_Test, coords_inside_image ) {
  vector<int> test;
  test = makeOverlapList( files_, inside_ );

  ASSERT_EQ(files_.size(), test.size()) << "Vectors are of unequal length.";
  for(  unsigned int i = 0; i < test.size(); ++i ) {
    EXPECT_EQ((int)i, test[i]) << "Coords did not overlap at index " << i;
  }
}

TEST_F( makeOverlapList_Test, coords_partiallly_inside_image ) {
  vector<int> test;
  test = makeOverlapList( files_, partial_ );

  ASSERT_EQ(files_.size(), test.size()) << "Vectors are of unequal length.";
  for(  unsigned int i = 0; i < test.size(); ++i ) {
    EXPECT_EQ((int)i, test[i]) << "Coords did not overlap at index " << i;
  }
}

TEST_F( makeOverlapList_Test, coords_outside_image ) {
  vector<int> test;
  test = makeOverlapList( files_, outside_ );

  ASSERT_TRUE( test.empty() ) << "Return vector should be empty.";
}

class makeOverlapListFromGeoTiff_Test : public ::testing::Test {
  protected:
  virtual void SetUp() {
    files_.push_back( "M111578606RE.10mpp.tif" );
    files_.push_back( "M111578606RE.10mpp.tif" );
    // MinimumLatitude      = 25.591191678781
    // MaximumLatitude      = 26.535020201925
    // MinimumLongitude     = 3.5042809131219
    // MaximumLongitude     = 3.6015317721328
    inside_(0) = 3.58;
    inside_(1) = 3.59;
    inside_(2) = 26;
    inside_(3) = 26.2;
    partial_(0) = 3.6;
    partial_(1) = 4.0;
    partial_(2) = 26;
    partial_(3) = 26.2;
    outside_(0) = 0;
    outside_(1) = 1;
    outside_(2) = 0;
    outside_(3) = 1;
  }

  virtual void TearDown() {}

  vector<string> files_;
  Vector4 inside_;
  Vector4 partial_;
  Vector4 outside_;
};

TEST_F( makeOverlapListFromGeoTiff_Test, coords_inside_image ) {
  vector<int> test;
  test = makeOverlapListFromGeoTiff( files_, inside_ );

  ASSERT_EQ(files_.size(), test.size()) << "Vectors are of unequal length.";
  for(  unsigned int i = 0; i < test.size(); ++i ) {
    EXPECT_EQ((int)i, test[i]) << "Coords did not overlap at index " << i;
  }
}

TEST_F( makeOverlapListFromGeoTiff_Test, coords_partiallly_inside_image ) {
  vector<int> test;
  test = makeOverlapListFromGeoTiff( files_, partial_ );

  ASSERT_EQ(files_.size(), test.size()) << "Vectors are of unequal length.";
  for(  unsigned int i = 0; i < test.size(); ++i ) {
    EXPECT_EQ((int)i, test[i]) << "Coords did not overlap at index " << i;
  }
}

TEST_F( makeOverlapListFromGeoTiff_Test, coords_outside_image ) {
  vector<int> test;
  test = makeOverlapListFromGeoTiff( files_, outside_ );

  ASSERT_TRUE( test.empty() ) << "Return vector should be empty.";
}

TEST( Save_and_ReadStatistics_Test, file_write_and_read ){
  fs::path p( "Save_and_ReadStatistics_Test.txt" );
  valarray<float> e(5);
  e[0] = 1.0;
  e[1] = 2.2;
  e[2] = 3.3;
  e[3] = 7.9;
  e[4] = 4.4;
  vector<float> b(4);
  b[0] = 0;
  b[1] = 2.5;
  b[2] = 5;
  b[3] = 7.5;
  vector<int> h(4);
  float minValid = 0;
  h[0] = static_cast<valarray<float> >(      e[ e>b[0] && e<=b[1] ] ).size();
  h[1] = static_cast<valarray<float> >(      e[ e>b[1] && e<=b[2] ] ).size();
  h[2] = static_cast<valarray<float> >(      e[ e>b[2] && e<=b[3] ] ).size();
  h[3] = static_cast<valarray<float> >(      e[ e>b[3] ]            ).size();
  int valid = static_cast<valarray<float> >( e[ e>minValid ]        ).size();

  SaveStatistics( p.string(), e, b );

  vector<int> test_h;
  float test_minE;
  float test_maxE;
  float test_avgE;
  int   test_valid;

  ReadStatistics( p.string(), test_h, &test_minE, &test_maxE, &test_avgE, &test_valid );

  ASSERT_EQ(h.size(), test_h.size()) << "Vectors are of unequal length.";
  for(  unsigned int i = 0; i < test_h.size(); ++i ) {
    EXPECT_EQ(h[i], test_h[i]) << "Histograms didn't have same values at index " << i 
                               << " truth: " << h[i] 
                               << " from file: " << test_h[i];
  }

  EXPECT_EQ( e.min(), test_minE )                   << "Minimum error differed.";
  EXPECT_EQ( e.max(), test_maxE )                   << "Maximum error differed.";
  EXPECT_NEAR( e.sum()/e.size(), test_avgE, 0.005 ) << "Average error differed.";
  EXPECT_EQ( valid, test_valid )                    << "Number of valid points differ.";

  fs::remove( p );
}

TEST( SaveOverlapList_Test, int_vector ){
  fs::path p("SaveOverlapList_Test.txt");
  vector<int> truth(4);
  truth[0] = 0;
  truth[1] = 1;
  truth[2] = 2;
  truth[3] = 3;

  SaveOverlapList( p.string(), truth );
  
  vector<int> test;
  test = ReadVectorFrom<int>( p );

  ASSERT_EQ(truth.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth.size(); ++i ) {
    EXPECT_EQ(truth[i], test[i]) << "Vectors truth and test differ at index " << i;
  } 

  fs::remove(p);
}

TEST( SaveOverlapList_Test, string_vector ){
  fs::path p("SaveOverlapList_Test.txt");
  vector<string> truth(4);
  truth[0] = "This";
  truth[1] = "is";
  truth[2] = "a";
  truth[3] = "test.";

  SaveOverlapList( p.string(), truth );
  
  vector<string> test;
  test = ReadVectorFrom<string>( p );

  ASSERT_EQ(truth.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth.size(); ++i ) {
    EXPECT_EQ(truth[i], test[i]) << "Vectors truth and test differ at index " << i;
  } 

  fs::remove(p);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
