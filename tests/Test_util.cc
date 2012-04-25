#include <boost/filesystem.hpp>

#include "../util.h"
#include "gtest/gtest.h"

namespace fs = boost::filesystem;

using namespace std;

TEST( util, ReadFileToVector ) {

  vector<string> truth;
  truth.push_back( "This" );
  truth.push_back( "is" );
  truth.push_back( "a" );
  truth.push_back( "test." );
  fs::path p ( "ReadFileToVector_temp.txt" );

  std::ofstream file( p.c_str() );
  ASSERT_TRUE( file ) << "Can't open file " << p ;
  
  for( vector<string>::const_iterator i = truth.begin(); 
       i != truth.end();
       ++i ){
    file << *i << std::endl;
  }
  file.close();

  vector<string> test;
  ReadFileList( p.string(), test );
  fs::remove( p );

  ASSERT_EQ(truth.size(), test.size()) << "Vectors are of unequal length";

  for( unsigned int i = 0; i < truth.size(); ++i ) {
    EXPECT_EQ(truth[i], test[i]) << "Vectors truth and test differ at index " << i;
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
