#include <boost/filesystem.hpp>

#include "../tracks.h"
#include "gtest/gtest.h"

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


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
