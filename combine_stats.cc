// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

//boost
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

//VW
#include <vw/Core.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>

//ATK
#include "util.h"

using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

int main( int argc, char *argv[] )
{
 int verbose;
 std::vector<std::string> inputStatsFiles;

 po::options_description general_options("Options");
 general_options.add_options()
    ("inputStatsFiles,s", po::value<std::vector<std::string> >(&inputStatsFiles))
    ("verbose,v", po::value<int>(&verbose)->default_value(1), 
		  "Verbosity level, zero emits no messages.")
    ("help,h", "Display this help message");
 

  po::options_description options("Allowed Options");
  options.add(general_options);
  
  po::positional_options_description p;
  p.add("inputStatsFiles", -1);
  
  std::ostringstream usage;
  usage << "Description: main code for the computation of the cummulative statistics" << std::endl << std::endl;
  usage << general_options << std::endl;
 
  po::variables_map vm;
  try
    {
      po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
      po::notify( vm );
    }
  catch ( po::error &e )
   {
     std::cout << "An error occured while parsing command line arguments.\n";
     std::cout << "\t" << e.what() << "\n\n";
     std::cout << usage.str();
     return 1;
   }
  
  if( vm.count("help") )
    {
      std::cerr << usage.str() << std::endl;
      return 1;
    }
  
  vector<string> statsFiles; 
  int numStatsInputFiles = inputStatsFiles.size();
 
  if(numStatsInputFiles < 1) {
    std::cerr << "Error: Must specify at least  one Stats file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }

  if( numStatsInputFiles == 1) {
     string inputFileExtension = GetFilenameExt(inputStatsFiles[0]);
     cout<<inputFileExtension<<endl;
     if (inputFileExtension.compare(string("txt")) ==0){
        ReadFileList(inputStatsFiles[0], statsFiles);
     }
     else{
       statsFiles.resize(1);
       statsFiles[0] = inputStatsFiles[0]; 
     }
     
  }
  if  ( numStatsInputFiles > 1) {
       statsFiles.resize(numStatsInputFiles);
       for (int i = 0; i < numStatsInputFiles; i++){
	 statsFiles[i] = inputStatsFiles[i]; 
         cout<<statsFiles[i]<<endl;
       }
  }
  
  int numStatsFiles = statsFiles.size();
  cout<<"numStatsFiles="<<numStatsFiles<<endl;

  vector<int> overallHist;
  overallHist.resize(5);
  float overallAvgError = 0;
  int overallNumValidPts = 0;
  float overallMinError = 100000000.0;
  float overallMaxError = -100000000.0;

  for (int i = 0; i < statsFiles.size(); i++){
    //TODO:read stats files
    float avgError;
    float minError;
    float maxError;
    int numValidPts;
    vector<int> hist;
    ReadStatistics (statsFiles[i], hist, &minError, &maxError, 
                                   &avgError, &numValidPts);
    /*
    for(int j = 0; j < 4; j++){
       cout<<hist[j]<<endl;
    }
    cout<<"minError="<<minError<<endl;
    cout<<"maxError="<<maxError<<endl;
    cout<<"avgError="<<avgError<<endl;
    cout<<"numValidPts="<<numValidPts<<endl;
    */

    overallAvgError = overallAvgError + avgError*numValidPts;
    overallNumValidPts = overallNumValidPts + numValidPts;
    for (int j = 0; j < 5; j++){
      overallHist[j] = overallHist[j]+hist[j];
    }
    if (minError<overallMinError){
      overallMinError = minError;
    }   
    if (maxError>overallMaxError){
      overallMaxError = maxError;
    }   
  }

  //print report
  overallAvgError = overallAvgError/overallNumValidPts;
  cout<<"overallAvgError="<<overallAvgError<<endl;
  cout<<"overallMinError="<<overallMinError<<endl;
  cout<<"overallMaxError="<<overallMaxError<<endl;
  cout<<"overallNumValidPts="<<overallNumValidPts<<endl;
  for (int i = 0; i < 5; i++){
    cout<<"overallHist["<<i<<"]="<<overallHist[i]<<endl;
  }

  return 0;
}


