// __BEGIN_LICENSE__
// Copyright (C) 2012 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>
#include "settings.h"

/*

This program is mostly a toy to exercise the Boost::program_options
routines, however, it also can generate a 'default' settings file,
based on the actual program defaults.  So when we do a release, we
should build and run this program to create an example settings
file to ship.

However, I don't think that including this program would be terribly
useful to end-users.

*/

using namespace std;
namespace po = boost::program_options;

int main( int argc, char* argv[] )
{
  po::options_description general_opts("Generic Options");
  general_opts.add_options()
               ("help,h", "Display this help message")
               ("settings-filename,s", 
                po::value<string>(),
                "File to read additional settings from.");

  po::options_description common_opts = getCommonSettings();

  po::options_description hidden_opts("Hidden options");
  hidden_opts.add_options()
         ("outfile", po::value<string>(), "File to write Default Settings to.");        
  po::positional_options_description p;
  p.add("outfile", -1);
  
  po::options_description all("Allowed options");
  all.add(general_opts).add(common_opts).add(hidden_opts);
  
  po::variables_map vm;
  po::store( po::command_line_parser(argc,argv).options(all).positional(p).run(), vm );

  if( vm.count("settings-filename") ){
    ifstream configFile( vm["settings-filename"].as<string>().c_str() );
    po::store( po::parse_config_file(configFile, common_opts), vm );
  }
  po::notify(vm);

  std::ostringstream usage;
  usage << "Description: default settings file writer and settings file reader." 
        << std::endl 
        << "    To write out a default settings file, foo.txt, do:"
        << std::endl
        << "        % ./settings_parser foo.txt"
        << std::endl << std::endl;
  usage << general_opts << std::endl;
  
  if( (argc < 2) or (vm.count("help")) ){
    std::cerr << usage.str() << std::endl;
    return 0;
  }

  if( vm.count("outfile") ){
    ofstream file( vm["outfile"].as<string>().c_str() );
    if( !file ) {
      cerr << "Can't open output file \"" << argv[0] << "\"" << endl;
      return 1;
      }
    file << CommonSettings_to_string( vm );
    file.close();
  }
  else { cout << CommonSettings_to_string( vm ); }
  return 0;
}
