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
//#include "weights.h"




int main( int argc, char *argv[] ) {

  //int N = 10;
  //simple_function_test(N);

  GlobalParams globalParams;
  //globalParams.reflectanceType = NO_REFL;
  globalParams.reflectanceType = LUNAR_LAMBERT;
  //globalParams.reflectanceType = LAMBERT;
  globalParams.slopeType = 1;
  globalParams.shadowThresh = 40;
  globalParams.albedoInitType = 1;
  globalParams.exposureInitType = 1;

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


  string inputCSVFilename; 
  string inputDEMFilename;
  string DRGFilename;  


  //inputCSVFilename = string("../data/Apollo15-LOLA/RDR_2E4E_25N27NPointPerRow_csv_table.csv"); 
  inputCSVFilename = string("../data/Apollo15-LOLA/1E8E_21N28N/RDR_1E8E_21N28NPointPerRow_csv_table.csv");  
  inputDEMFilename = string("../data/Apollo15-DEM/1134_1135-DEM.tif");

  DRGFilename = string("../data/Apollo15-DRG/1134_1135-DRG.tif");  
  //DRGFilename = string("../data/Apollo15-DRG/AS15-M-1134_map.tif");  


  string DRGFilenameNoPath = sufix_from_filename(DRGFilename);

  string imgPtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "img.txt";    
  string reflectancePtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "refl.txt";  
  string syntImgPtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "synth.txt";  
  string altitudePtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "alt.txt";
  string demPtsFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "dem.txt";
  string lolaTracksFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "_lola.tif";  
  string outFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "_results.tif";  

  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);
  int analyseFlag = 0;//1; 
  int maxNumIter = 10;
  int maxNumStarts = 160;
  vector<Vector<float, 6> >init_d_array;
  init_d_array.resize(maxNumStarts);
  for (int i = 0; i < maxNumStarts; i++){
    init_d_array[i][0] = 1.0;
    init_d_array[i][1] = 0.0;
    init_d_array[i][2] = (i-maxNumStarts/2)*5;
    init_d_array[i][3] = 0.0;
    init_d_array[i][4] = 1.0;
    init_d_array[i][5] = 0.0;//(i-maxNumStarts/2)*25;
  }    

  vector<int> trackIndices;
  trackIndices.resize(trackPts.size());
  for (int i = 0; i < trackPts.size(); i++){
    trackIndices[i] = i;
  }

  vector<Vector<float, 6> > final_d_array;
  vector<float> error_array;
  error_array.resize(maxNumStarts);
  final_d_array.resize(maxNumStarts);


  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);


  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  //get the true image points
  cout << "GetAllPtsFromImage..." << endl; 
  GetAllPtsFromImage(trackPts, interpDRG, DRGGeo);

  cout << "ComputeTrackReflectance..." << endl;
  ComputeAllReflectance(trackPts, modelParams, globalParams);
  float scaleFactor = ComputeScaleFactor(trackPts);

  if (analyseFlag == 1){

    SaveImagePoints(trackPts, 3, imgPtsFilename);
    SaveAltitudePoints(trackPts, 3, altitudePtsFilename);
    SaveReflectancePoints(trackPts, 1.0, reflectancePtsFilename);
    SaveReflectancePoints(trackPts, scaleFactor, syntImgPtsFilename);
    SaveDEMPoints(trackPts, inputDEMFilename, demPtsFilename);

    int numVerPts = 6000;
    int numHorPts = 6000;
    MakeGrid(trackPts, numVerPts, numHorPts, lolaTracksFilename, trackIndices);
  }

  int edge = 10;
  float take_p = 0.10;
  int num_valid = 0;
  float take_thresh = 0.0;
  string s_weight_name = "../results/weights_corregister_prd.txt ";
  cout << "Calling weight_track_pts... "<< endl;
  weight_track_pts( trackPts, edge, take_p, num_valid, take_thresh, s_weight_name );

  printf("Weight calc: edge = %d, take_p = %f, num_valid = %d, take_thresh = %d\n\n", edge, take_p, num_valid, take_thresh);


  //return matching error and transform
  cout << "Calling UpdateMatchingParams ..."<< endl;
  UpdateMatchingParams(trackPts, DRGFilename, 
      modelParams, globalParams,maxNumIter,  
      init_d_array, final_d_array, error_array);

  int bestResult = 0;
  float smallestError = error_array[0];
  for (int index = 0; index < init_d_array.size(); index++){
    printf("OUT %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
        index, error_array[index], 
        final_d_array[index](0), final_d_array[index](1), 
        final_d_array[index](2), final_d_array[index](3),
        final_d_array[index](4), final_d_array[index](5));
    if  (error_array[index] < smallestError){
      smallestError = error_array[index]; 
      bestResult = index;
    }      
  }    

  //write results to image outside matching

  printf("bestResult = %d\n", bestResult);
  ShowFinalTrackPtsOnImage(trackPts, final_d_array[bestResult], 
      trackIndices, DRGFilename, outFilename);

  cout << "UpdateMatchingParams finshed." << endl;

  return 0;
}


















