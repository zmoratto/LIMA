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
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>

//ATK
#include "../util.h"

using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

int main( int argc, char *argv[] )
{
 int verbose;

 std::vector<std::string> inputFiles;
 std::string resDir;
 float north, south, east, west;

 po::options_description general_options("Options");
 general_options.add_options()
    ("inputFiles,d", po::value<std::vector<std::string> >(&inputFiles))
    ("north,n", po::value<float> (&north))
    ("south,s", po::value<float> (&south))
    ("east,e", po::value<float> (&east)) 
    ("west,w", po::value<float> (&west))
    ("results-directory,r", 
		po::value<std::string>(&resDir)->default_value("../results"), 
		"results directory.") // Currently a no-op until we implement it below.
    ("verbose,v", 
		po::value<int>(&verbose)->default_value(1), 
		"Verbosity level, zero emits no messages.")
    ("help,h", "Display this help message");
  

  po::options_description hidden_options("");

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("inputFiles", -1);

  std::ostringstream usage;
  usage << "Description: main application for georeferenced image selection within a lonlat bounding box" << std::endl << std::endl;
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

  vector<string> imageFiles; 
  int numInputFiles = inputFiles.size();

  if(numInputFiles < 1) {
    std::cerr << "Error: Must specify at least  one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }

  //determine if inputFiles is a text file containing a list of image files, one image file or a set of image files
  //by reading the file extension. A text file containing the image list *must* have extension .txt 
  imageFiles = AccessDataFilesFromInput(inputFiles);
 
  //Set up VW logging
  vw_log().console_log().rule_set().add_rule(vw::InfoMessage,"*");

  int numImageFiles = imageFiles.size();
  cout<<"numImageFiles="<<numImageFiles<<endl;
 
  cout<<"Selecting the overlapping DEMs ..."<<endl;
  std::vector<int> overlapIndices;
  string overlapListFilename = resDir+string("/")+string("list.txt");

  Vector4 lon_lat_bb;
  lon_lat_bb[0]=west;  //minlon
  lon_lat_bb[1]=east;  //maxlon
  lon_lat_bb[2]=south; //minlat
  lon_lat_bb[3]=north; //maxlat
  overlapIndices = makeOverlapListFromGeoTiff(imageFiles, lon_lat_bb);
  //SaveOverlapList(overlapListFilename, overlapIndices);
  
  std::vector<std::string> overlapFilenames;
  for (int i= 0; i < overlapIndices.size(); i++){
    overlapFilenames.push_back(imageFiles[overlapIndices[i]]);
  }
  SaveOverlapList(overlapListFilename, overlapFilenames);
  
  //PrintOverlapList(overlapIndices);
  for (int i= 0; i < overlapIndices.size(); i++){
    cout<<imageFiles[overlapIndices[i]]<<endl;
  }
  cout <<"overlapListFilenamae="<<overlapListFilename<<endl;

  return 0;
}

