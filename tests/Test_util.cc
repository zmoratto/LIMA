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

TEST( util, ReadFileToVector ) {
  vector<string> truth(4);
  truth[0] = "This";
  truth[1] = "is";
  truth[2] = "a";
  truth[3] = "test.";
  fs::path p ( "ReadFileToVector_temp.txt" );

  std::ofstream file( p.c_str() );
  ASSERT_TRUE( file ) << "Can't open file " << p ;
  test_write_vector_to_stream( truth, file );
  file.close();

  vector<string> test;
  ReadFileList( p.string(), test );
  fs::remove( p );

  ASSERT_EQ(truth.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth.size(); ++i ) {
    EXPECT_EQ(truth[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST( util, AccessDataFilesFromInput ) {
  vector<string> truth1(1,"File1");
  vector<string> truth2(3);
  truth2[0] = "File01"; 
  truth2[1] = "File02"; 
  truth2[2] = "File03.data"; 

  fs::path truth1_p( "AccessDataFilesFromInput_truth1.txt" );
  fs::path truth2_p( "AccessDataFilesFromInput_truth2.txt" );

  vector<string> list1(1, truth1_p.string() );
  // vector<string> list2(4);
  // list2[0] = "File_A"; 
  // list2[1] = "File_B"; 
  // list2[2] = "File_C.data"; 
  // list2[3] = truth2_p.string(); 
  // vector<string> truth3( list2.begin(), list2.end()-1 );
  // truth3.insert( truth3.end(), truth2.begin(), truth2.end() );

  // Tests for plain vectors
  vector<string> test;
  test = AccessDataFilesFromInput( truth1 );
  ASSERT_EQ(truth1.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth1.size(); ++i ) {
    EXPECT_EQ(truth1[i], test[i]) << "Vectors truth and test differ at index " << i;
  }

  test.clear();
  test = AccessDataFilesFromInput( truth2 );
  ASSERT_EQ(truth2.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth2.size(); ++i ) {
    EXPECT_EQ(truth2[i], test[i]) << "Vectors truth and test differ at index " << i;
  }

  // Write out some files to test the file reading part
  std::ofstream file0( truth1_p.c_str() );
  ASSERT_TRUE( file0 ) << "Can't open file " << truth1_p ;
  test_write_vector_to_stream( truth1, file0 );
  file0.close();

  std::ofstream file1( truth2_p.c_str() );
  ASSERT_TRUE( file1 ) << "Can't open file " << truth2_p ;
  test_write_vector_to_stream( truth2, file1 );
  file1.close();

  // Test for vectors that contain file list names
  test.clear();
  test = AccessDataFilesFromInput( list1 );
  ASSERT_EQ(truth1.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth1.size(); ++i ) {
    EXPECT_EQ(truth1[i], test[i]) << "Vectors truth and test differ at index " << i;
  }

  /* Run this test when input files can contain .txt files that can load in other lists.
  test.clear();
  test = AccessDataFilesFromInput( list2 );
  ASSERT_EQ(truth3.size(), test.size()) << "Vectors are of unequal length";
  for( unsigned int i = 0; i < truth3.size(); ++i ) {
    EXPECT_EQ(truth3[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
  */

  //Clean up files
  fs::remove( truth1_p );
  fs::remove( truth2_p );
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
