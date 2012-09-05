// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef LIDAR2DEM_CONFIG_H
#define LIDAR2DEM_CONFIG_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/lexical_cast.hpp>

#include <vw/Core.h>
#include <vw/Math.h>

using namespace vw;
using namespace vw::math;

using namespace std;
/*
inline std::ostream& operator<< ( std::ostream& os, const CoregistrationParams& cp )
	{
	os	<< "  MATCHING_MODE: " 
		<< boost::lexical_cast<std::string>(cp.matchingMode)			<< endl
		<< "  REFLECTANCE_TYPE: " 
		<< boost::lexical_cast<std::string>(cp.reflectanceType)			<< endl
		<< "  ANALYSE_FLAG: " 
		<< boost::lexical_cast<std::string>(cp.analyseFlag)				<< endl
		<< "  USE_REFLECTANCE_FEATURES: " 
		<< boost::lexical_cast<std::string>(cp.useReflectanceFeatures)	<< endl
		<< "  TOP_PERCENT_FEATURES: " 
		<< boost::lexical_cast<std::string>(cp.topPercentFeatures)		<< endl
		<< "  samplingStep: " << cp.samplingStep << endl
		<< "  matchWindowHalfSize: "<< cp.matchWindowHalfSize			<< endl
		<< "  MAX_NUM_ITER: " 
		<< boost::lexical_cast<std::string>(cp.maxNumIter)				<< endl
		<< "  MAX_NUM_STARTS: " 
		<< boost::lexical_cast<std::string>(cp.maxNumStarts)			<< endl
	  //<< "  DISPLAY_RESULTS: " 
	  //<< boost::lexical_cast<std::string>(cp.displayResults)			<< endl
		<< "  NO_DATA_VAL: " 
		<< boost::lexical_cast<std::string>(cp.noDataVal)				<< endl
		<< "  MIN_CONV_THRESH: " 
		<< boost::lexical_cast<std::string>(cp.minConvThresh)			<< endl;
		
	return os;
	};
*/
bool ReadConfigFile(string config_filename, struct CoregistrationParams *settings);

void SaveVectorToFile(vector<float> v, string filename);
#endif /* LIDAR2DEM_CONFIG_H */
