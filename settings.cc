// __BEGIN_LICENSE__
// Copyright (C) 2012 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <fstream>
#include <iostream>
#include <boost/any.hpp>
#include <boost/program_options.hpp>
#include "settings.h"

/*

The intent with these bits of code is to create a core set of shared
settings that all the tools in AlignTK would use, and a basic settings
file.  Individual tools will want to define additional settings specific
to their needs, as well.

*/

using namespace std;
namespace po = boost::program_options;

po::options_description getCommonSettings(){
  po::options_description desc;
  // If you add or change an option here, make sure to add it below in CommonSettings_to_string()
  desc.add_options() 
      ("MATCHING_MODE", 
       po::value<int>()->default_value(1), 
       "0 - no matching, 1 - affine 2D, 2 - ICP 3D")
      ("SAMPLING_STEP_H",
       po::value<int>()->default_value(8),
       "The horizontal step size through the DEM or image data.")
      ("SAMPLING_STEP_V",
       po::value<int>()->default_value(8),
       "The vertical step size through the DEM or image data.")
      ("MATCH_WINDOW_H",
       po::value<int>()->default_value(5),
       "The half width of the matching window.")
      ("MATCH_WINDOW_V",
       po::value<int>()->default_value(5),
       "The half height of the matching window.")
      ("MAX_NUM_ITER",
       po::value<int>()->default_value(3),
       "The maximum number of iterations the program will make.")
      ("MAX_NUM_STARTS",
       po::value<int>()->default_value(160),
       "??")
      ("NO_DATA_VAL",
       po::value<double>()->default_value(-10000),
       "No data value for images or DEMs that don’t have a no data speciﬁed in their format.")
      ("CONV_THRESH",
       po::value<float>()->default_value(0.01),
       "The threshold under which absolute value at consecutive iterations determine the algorithm convergence.");

  return desc;
}

string CommonSettings_to_string( const po::variables_map& vm ){
 /* I initially thought that just iterating through vm and printing as I went 
 * would be more elegant [see VariablesMap_to_string() below].  However, it just prints
 * everything in the vm, and there's no control over the order [although it does seem
 * to be alphabetical, I don't know if that's guaranteed].  So I just went 
 * with this manual approach.
 */
  ostringstream s;
  s << "MATCHING_MODE = "   << vm["MATCHING_MODE"].as<int>()   << endl;
  s << "SAMPLING_STEP_H = " << vm["SAMPLING_STEP_H"].as<int>() << endl;
  s << "SAMPLING_STEP_V = " << vm["SAMPLING_STEP_V"].as<int>() << endl;
  s << "MATCH_WINDOW_H = "  << vm["MATCH_WINDOW_H"].as<int>()  << endl;
  s << "MATCH_WINDOW_V = "  << vm["MATCH_WINDOW_V"].as<int>()  << endl;
  s << "MAX_NUM_ITER = "    << vm["MAX_NUM_ITER"].as<int>()    << endl;
  s << "MAX_NUM_STARTS = "  << vm["MAX_NUM_STARTS"].as<int>()  << endl;
  s << "NO_DATA_VAL = "     << vm["NO_DATA_VAL"].as<double>()  << endl;
  s << "CONV_THRESH = "     << vm["CONV_THRESH"].as<float>()   << endl;
  return s.str();
}


string VariablesMap_to_string( const po::variables_map& vm ){
  ostringstream s;
  map<std::string, po::variable_value>::const_iterator iter;
  for( iter = vm.begin(); iter != vm.end(); ++iter ){
    if( typeid(int) == iter->second.value().type() ){
      s << iter->first << " = " << boost::any_cast<int>(iter->second.value()) << endl;
    }
    else if( typeid(float) == iter->second.value().type() ){
      s << iter->first << " = " << boost::any_cast<float>(iter->second.value()) << endl;
    }
    else if( typeid(double) == iter->second.value().type() ){
      s << iter->first << " = " << boost::any_cast<double>(iter->second.value()) << endl;
    }
    else if( typeid(string) == iter->second.value().type() ){
      s << iter->first << " = " << boost::any_cast<string>(iter->second.value()) << endl;
    }
    else {
      s << "# The value of " << iter->first << " couldn't be auto-converted." << endl;
    }
  }
  return s.str();
}
