// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/lexical_cast.hpp>
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Photometry.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>


using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;

#include <math.h>
#include "coregister.h"


bool ReadConfigFile(string config_filename, struct CoregistrationParams *settings)
{
  int MAX_LENGTH = 5000;
  char line[MAX_LENGTH];
  ifstream configFile(config_filename.c_str());

  if (configFile.is_open()){
    
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MATCHING_MODE %d\n", &(settings->matchingMode));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "REFLECTANCE_TYPE %d\n", &(settings->reflectanceType));
   
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "ANALYSE_FLAG %d\n", &(settings->analyseFlag));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "USE_REFLECTANCE_FEATURES %d\n", &(settings->useReflectanceFeatures));
   
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "TOP_PERCENT_FEATURES %d\n", &(settings->topPercentFeatures));

    configFile.getline(line, MAX_LENGTH);
    int samplingStepX, samplingStepY;
    sscanf(line, "SAMPLING_STEP %d %d\n", &samplingStepX, &samplingStepY);
    settings->samplingStep(0) = samplingStepX;
    settings->samplingStep(1) = samplingStepY; 

    configFile.getline(line, MAX_LENGTH);
    int windowSizeX, windowSizeY;
    sscanf(line, "MATCH_WINDOW %d %d\n", &windowSizeX, &windowSizeY);
    settings->matchWindowHalfSize(0) = windowSizeX;
    settings->matchWindowHalfSize(1) = windowSizeY;
   
 
    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MAX_NUM_ITER %d\n", &(settings->maxNumIter));

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "MAX_NUM_STARTS %d\n", &(settings->maxNumStarts)); 

    //configFile.getline(line, MAX_LENGTH);
    //sscanf(line, "DISPLAY_RESULTS %d\n", &(settings->displayResults)); 

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "NO_DATA_VAL %lf\n", &(settings->noDataVal)); 

    configFile.getline(line, MAX_LENGTH);
    sscanf(line, "CONV_THRESH %f\n", &(settings->minConvThresh)); 
    
    configFile.close();

    return true;
  }
  else{
    settings->matchingMode = 1; // LIMA;//LIDEM
    settings->reflectanceType = 3; //LUNAR_LAMBERT;// NO_REFL;//LAMBERT;
    settings->analyseFlag = 0;
    settings->useReflectanceFeatures = 0;
    settings->topPercentFeatures = 10;
    settings->samplingStep(0) = 8;
    settings->samplingStep(1) = 8;
    settings->matchWindowHalfSize(0) = 5;
    settings->matchWindowHalfSize(1) = 5;
    settings->maxNumIter = 3;
    settings->maxNumStarts = 160;
    //settings->displayResults = 0; 
    settings->noDataVal = -10000;
    settings->minConvThresh = 0.01;

	return false;
  }
}

/*
//this function will be removed
int ReadModelParamsFile(string modelParamsFilename, struct ModelParams *params)
{

 int MAX_LENGTH = 5000;
 char line[MAX_LENGTH];
 ifstream modelParamsFile((char*)modelParamsFilename.c_str());

  if (modelParamsFile.is_open()){
    printf("MODEL PARAMS FILE FOUND\n");
    
    modelParamsFile.getline(line, MAX_LENGTH);

    modelParamsFile.getline(line, MAX_LENGTH);
    sscanf(line, "LIGHT_POS %lf %lf %lf\n", &(params->sunPosition[0]), &(params->sunPosition[1]), &(params->sunPosition[2])  );

    modelParamsFile.getline(line, MAX_LENGTH);
    sscanf(line, "VIEWER_POS %lf %lf %lf\n", &(params->spacecraftPosition[0]), &(params->spacecraftPosition[1]), &(params->spacecraftPosition[2]));

    params->exposureTime = 1.0;
    //params->rescalingParams[0] = 1;
    //params->rescalingParams[1] = 0;

    params->sunPosition[0] = 1000*(params->sunPosition[0]);
    params->sunPosition[1] = 1000*(params->sunPosition[1]);
    params->sunPosition[2] = 1000*(params->sunPosition[2]);

    params->spacecraftPosition[0] = 1000*(params->spacecraftPosition[0]);
    params->spacecraftPosition[1] = 1000*(params->spacecraftPosition[1]);
    params->spacecraftPosition[2] = 1000*(params->spacecraftPosition[2]);
    return (0);
  }

  return(1);
}
*/

void PrintModelParams(struct ModelParams *params)
{
  cout<<"LIGHT_POS[0]"<<params->sunPosition[0]<<endl;
  cout<<"LIGHT_POS[1]"<<params->sunPosition[1]<<endl;
  cout<<"LIGHT_POS[2]"<<params->sunPosition[2]<<endl;

  cout<<"VIEWER_POS[0]"<<params->spacecraftPosition[0]<<endl;
  cout<<"VIEWER_POS[1]"<<params->spacecraftPosition[1]<<endl;
  cout<<"VIEWER_POS[2]"<<params->spacecraftPosition[2]<<endl;
}



void SaveVectorToFile(vector<float> v, string filename)
{
 

  FILE *fp;

  fp = fopen(filename.c_str(), "w");

  for (unsigned int i = 0; i < v.size()-1; i++){
    fprintf(fp, "%f\n", v[i]);
  }
  fprintf(fp, "%f", v[v.size()-1]);
  fclose(fp);
}



















