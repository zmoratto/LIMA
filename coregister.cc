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
#include "tracks.h"
#include "io.h"
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"


int main( int argc, char *argv[] ) {


  string inputCSVFilename; 
  string inputDEMFilename;
  string inputDRGFilename;  

  if (argc == 5){ 

     if (strcmp(argv[1], "-i") == 0){
         inputDRGFilename = string(argv[2]);
     }
     else{
         if (strcmp(argv[1], "-d") == 0){
             inputDEMFilename = string(argv[2]);
         }
         else{
	    cout<<"Usage: coregister -i DRGFilename -d DEMFilename -l CSVFilename"<<endl;
	    return -1;
         }
     }

     if (strcmp(argv[3], "-l") == 0){
        inputCSVFilename = string(argv[4]);
     }
     else{
          cout<<"Usage: coregister -i DRGFilename -d DEMFilename -l CSVFilename"<<endl;
          return -1;
     }
 
  }
  else{
    cout<<"Usage: coregister -i DRGFilename -d DEMFilename -l CSVFilename"<<endl;
    return -1;
  }

  struct CoregistrationParams settings;
  string configFilename="coregister_settings.txt";
  ReadConfigFile((char*)configFilename.c_str(), &settings);
  PrintGlobalParams(&settings);

  //these values will be read from files
  ModelParams modelParams;
  modelParams.exposureTime = 1.0;
  modelParams.rescalingParams[0] = 1;
  modelParams.rescalingParams[1] = 0;

  modelParams.sunPosition[0] = 1000*72987509.682619;//*sunPositions[i][0];
  modelParams.sunPosition[1] = 1000*133319340.07726;//*sunPositions[i][1];
  modelParams.sunPosition[2] = 1000*555820.93952321;//*sunPositions[i][2];

  modelParams.spacecraftPosition[0] = 1000*1668.1656423675;
  modelParams.spacecraftPosition[1] = 1000*127.53606522819;
  modelParams.spacecraftPosition[2] = 1000*774.17340580747;

  
  //create the results directory and prepare the output filenames - START
  system("mkdir ../results");

  string DRGFilenameNoPath = sufix_from_filename(inputDRGFilename);

  string imgPtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "img.txt";    
  string reflectancePtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "refl.txt";  
  string syntImgPtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "synth.txt";  
  string altitudePtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "alt.txt";
  string demPtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "dem.txt";
  string lolaTracksFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "_lola.tif";  
  string outFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "_results.tif";  
  string lolaFeaturesFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "_features_lola.txt";  
  string matchResultsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "match_results"  +".txt";  

  //create the results directory and prepare the output filenames - END

  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);

  vector<Vector<float, 6> >initTransfArray;
  initTransfArray.resize(settings.maxNumStarts);
  for (int i = 0; i < settings.maxNumStarts; i++){
    initTransfArray[i][0] = 1.0;
    initTransfArray[i][1] = 0.0;
    initTransfArray[i][2] = (i-settings.maxNumStarts/2)*5;
    initTransfArray[i][3] = 0.0;
    initTransfArray[i][4] = 1.0;
    initTransfArray[i][5] = 0.0;//(i-maxNumStarts/2)*25;
  }    

  vector<int> trackIndices;
  trackIndices.resize(trackPts.size());
  for (int i = 0; i < trackPts.size(); i++){
    trackIndices[i] = i;
  }

  vector<Vector<float, 6> > finalTransfArray;
  vector<float> errorArray;
  errorArray.resize(settings.maxNumStarts);
  finalTransfArray.resize(settings.maxNumStarts);


  //initialization step for LIMA - START  
  DiskImageView<PixelGray<uint8> >   DRG(inputDRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, inputDRGFilename);


  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
									ConstantEdgeExtension()),
							    BilinearInterpolation());
  
  //get the true image points
  cout << "GetAllPtsFromImage..." << endl; 
  GetAllPtsFromImage(trackPts, interpDRG, DRGGeo);
  
  cout << "ComputeTrackReflectance..." << endl;
  ComputeAllReflectance(trackPts, modelParams, settings);
  //initialization step for LIMA - END  

  //TO DO: initialization step for LIDEM

  if (settings.analyseFlag == 1){
    float scaleFactor = ComputeScaleFactor(trackPts);
    SaveImagePoints(trackPts, 3, imgPtsFilename);
    SaveAltitudePoints(trackPts, 3, altitudePtsFilename);
    SaveReflectancePoints(trackPts, 1.0, reflectancePtsFilename);
    SaveReflectancePoints(trackPts, scaleFactor, syntImgPtsFilename);
    SaveDEMPoints(trackPts, inputDEMFilename, demPtsFilename);

    int numVerPts = 6000;
    int numHorPts = 6000;
    MakeGrid(trackPts, numVerPts, numHorPts, lolaTracksFilename, trackIndices);
  }

 
  if (settings.useLOLAFeatures){

    cout << "Computing the LOLA features and weights ... ";
    int halfWindow = 10;
    float topPercent = 0.10;
    ComputeWeights( trackPts, halfWindow, topPercent, lolaFeaturesFilename);
    cout<<"done."<<endl;

  }

  if (settings.matchingMode == LIMA){
    //return matching error and transform
    cout << "UpdateMatchingParams ..."<< endl;
    UpdateMatchingParamsLIMA_MP(trackPts, inputDRGFilename, 
				modelParams, settings,  
				initTransfArray, finalTransfArray, errorArray);
  }

  if (settings.matchingMode == LIDEM){
    //return matching error and transform
    cout << "UpdateMatchingParams ..."<< endl;
    UpdateMatchingParamsLIDEM_MP(trackPts, inputDEMFilename, 
			         modelParams, settings,  
			         initTransfArray, finalTransfArray, errorArray);
  }

  int bestResult = 0;
  float smallestError = errorArray[0];
  for (int index = 0; index < initTransfArray.size(); index++){
    printf("OUT %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
	   index, errorArray[index], 
	   finalTransfArray[index](0), finalTransfArray[index](1), 
	   finalTransfArray[index](2), finalTransfArray[index](3),
	   finalTransfArray[index](4), finalTransfArray[index](5));
    if  (errorArray[index] < smallestError){
      smallestError = errorArray[index]; 
      bestResult = index;
    }      
  }    
  cout<<"bestResult= "<<bestResult<<endl;


  //write finalTransfArray and errorArray to file
  SaveMatchResults(finalTransfArray, errorArray, matchResultsFilename);

  if (settings.displayResults){
    //write results to image outside matching
    ShowFinalTrackPtsOnImage(trackPts, finalTransfArray[bestResult], 
                             trackIndices, inputDRGFilename, outFilename);
  }
  cout << "UpdateMatchingParams done." << endl;

  return 0;
}


















