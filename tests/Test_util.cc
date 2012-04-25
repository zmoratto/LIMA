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

    p_ = "ReadFileToVector_temp.txt";

    std::ofstream file( p_.c_str() );
    ASSERT_TRUE( file ) << "Can't open file " << p_ ;
    test_write_vector_to_stream( truth_, file );
    file.close();
  }

  virtual void TearDown() {
    fs::remove( p_ );
  }

  vector<string> truth_;
  fs::path p_;
};

TEST_F( ReadFileToVectorTest, returnsString ) {
  vector<string> test = ReadFileList( p_.string() );

  ASSERT_EQ(truth_.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth_.size(); ++i ) {
    EXPECT_EQ(truth_[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

TEST_F( ReadFileToVectorTest, returnsVoid ) {
  vector<string> test;
  ReadFileList( p_.string(), test );

  ASSERT_EQ(truth_.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth_.size(); ++i ) {
    EXPECT_EQ(truth_[i], test[i]) << "Vectors truth and test differ at index " << i;
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


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
