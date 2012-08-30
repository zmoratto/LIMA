// __BEGIN_LICENSE__
// Copyright (C) 2012 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include "tracks.h"
#include "util.h"

using namespace std;

int main( int argc, char* argv[] )
{
po::options_description opts("Program Options");
opts.add_options()
  ("help,h", "Display this help message")
  ("lidar-filename,l", po::value<std::string>())
  ("output-prefix,o", po::value<std::string>())
  ("output-directory,r", po::value<std::string>()->default_value("./"));

po::positional_options_description p;
p.add("lidar-filename", -1);

std::ostringstream usage;
  usage << "Description: Takes a LOLA *PointPerRow_csv_table.csv file and splits it into several files, one per track." << std::endl << std::endl;
  usage << opts << std::endl;

po::variables_map vm;
try{
  po::store( po::command_line_parser(argc,argv).options(opts).positional(p).run(), vm );
}
catch( po::error& e ){
  cout << "An error occured while parsing command line arguments." << endl;
  cout << "\t" << e.what() << endl << endl;
  cout << usage.str();
  return 1;
}

if( vm.count("help") ){
  std::cerr << usage.str() << endl;
  return 0;
}

if( !vm.count("lidar-filename") ){
  cerr << "Error: Must specify a Lidar file!" << endl << endl;
  cerr << usage.str();
  return 1;
}

fs::path tracks_path( vm["lidar-filename"].as<string>() );
fs::path output_path( vm["output-directory"].as<string>() );
string output_prefix;
if( vm.count( "output-prefix" ) ){
  output_prefix = vm["output-prefix"].as<string>();
}
else {
  output_prefix = tracks_path.stem().string();
  output_prefix.append("_");
}

// if output_path doesn't exist, make it.
if( vm["output-directory"].as<string>() != "./" ){
  if( fs::exists(output_path) ){
    if( !fs::is_directory(output_path) ){
      cerr << output_path << " is not a directory." << endl;
      return 1;
    }
  }
  else{
    if( fs::create_directory( output_path ) == false ){
      cerr << "Could not create the " << output_path << " directory." << endl;
      return 1;
    }
  }
}

// Read in LOLA file
// I'd like to just use LOLAFileRead() here, but it would require retooling 
// the method to retain the full line which was read in for each shot, which 
// seems excessive.  Instead, I just reimplement the minimum here.  So if the
// logic changes there, it will probably need to change here.

ifstream file( tracks_path.c_str() );
  if( !file ) {
    cerr << "Unable to open track file " << tracks_path << endl;
    return 1;
  }

boost::char_separator<char> sep(" ,");
pointCloud currPt;
pointCloud prevPt;
string line;
ofstream outfile;
while( getline(file, line, '\n') ) {
  boost::tokenizer<boost::char_separator<char> > tokens( line, sep );
  boost::tokenizer<boost::char_separator<char> >::iterator token = tokens.begin();
  istringstream date_s( *token );
  date_s >> setw(4) >> currPt.year
	     >> ignoreOne // -
	     >> setw(2) >> currPt.month
	     >> ignoreOne // -
	     >> setw(2) >> currPt.day
	     >> ignoreOne // T
	     >> setw(2) >> currPt.hour
	     >> ignoreOne // :
	     >> setw(2) >> currPt.min
	     >> ignoreOne // :
	     >> currPt.sec;
  if( currPt.year <= 0 ) { continue; }

  if( prevPt.year == 0 || isTimeDiff(prevPt, currPt, 3000) ) { //new track
    if( outfile.is_open() ){
      outfile.close();
    }
    string outname( output_prefix );
    outname.append( boost::algorithm::replace_all_copy(currPt.time(),":","-") );
    outname.append( ".csv" );
    fs::path p( output_path/outname );
    outfile.open( p.c_str() );
  }
  outfile << line << endl;

  prevPt = currPt;
} 

file.close();

return 0;
}
