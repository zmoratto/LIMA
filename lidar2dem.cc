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
#include "icp.h"

int main( int argc, char *argv[] ) {


  string inputCSVFilename; 
  string inputDEMFilename;


  if (argc == 5){ 

    
     if (strcmp(argv[1], "-d") == 0){
         inputDEMFilename = string(argv[2]);
     }
     else{
	    cout<<"Usage: coregister -d DEMFilename -l CSVFilename"<<endl;
	    return -1;
     }
    
     if (strcmp(argv[3], "-l") == 0){
        inputCSVFilename = string(argv[4]);
     }
     else{
          cout<<"Usage: coregister -d DEMFilename -l CSVFilename"<<endl;
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

   
  ModelParams modelParams;
  string modelParamsFilename="light_viewer_pos.txt";
  ReadModelParamsFile(modelParamsFilename, &modelParams);
  PrintModelParams(&modelParams);
  
  //create the results directory and prepare the output filenames - START
  system("mkdir ../results");

  string DEMFilenameNoPath = sufix_from_filename(inputDEMFilename);

  string altitudePtsFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "alt.txt";
  string demPtsFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "dem.txt";
  string lolaTracksFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "_lola.tif";  
  string outFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "_results.tif";  
  string lolaFeaturesFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "_features_lola.txt";  
  string matchResultsFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "match_results"  +".txt";  

  //create the results directory and prepare the output filenames - END

  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);
 
  vector<Vector<float, 6> > finalTransfArray;
  vector<float> errorArray;
  errorArray.resize(settings.maxNumStarts);
  finalTransfArray.resize(settings.maxNumStarts);


  //initialization step for LIDE - START
  DiskImageView<PixelGray<float> >   DEM(inputDEMFilename);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, inputDEMFilename);


  ImageViewRef<PixelGray<float> >   interpDEM = interpolate(edge_extend(DEM.impl(),
									ConstantEdgeExtension()),
							    BilinearInterpolation());


  GetAllPtsFromDEM(trackPts, interpDEM, DEMGeo);
  //initialization step for LIDEM - END
  /*
  if (settings.analyseFlag == 1){
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
  */
  Vector3 currTranslation;
  Matrix<float, 3,3 > currRotation;
  vector<Vector3> featureArray;//lidarData
  vector<Vector3> modelArray;//DEM
  
  
  //copy info to featureArray and modelArray
  for(int k = 0; k < trackPts.size();k++){
     for(int i = 0; i < trackPts[k].size(); i=i+100){
       
       if ((trackPts[k][i].valid == 1) && (trackPts[k][i].DEMPt[2].valid==1)){
	 Vector3 model;
	 Vector3 feature;
         
         //this is the LOLA data
	 model[0] = trackPts[k][i].LOLAPt[2].coords(0); 
	 model[1] = trackPts[k][i].LOLAPt[2].coords(1);   
	 model[2] = trackPts[k][i].LOLAPt[2].coords(2);

         if ((model[0] >1e-100) && (model[1] >1e-100) && (model[2]>1e-100)){
	   feature[0] = trackPts[k][i].LOLAPt[2].coords(0); 
	   feature[1] = trackPts[k][i].LOLAPt[2].coords(1); 
	   feature[2] = trackPts[k][i].DEMPt[2].val;

	   //copy the model and features into cartesian coordinates
	   feature = DEMGeo.datum().geodetic_to_cartesian(feature);
	   model = DEMGeo.datum().geodetic_to_cartesian(model);
         	 
	   featureArray.push_back(feature);
	   modelArray.push_back(model);
	 }
       }
     }
  }

  cout<<modelArray.size()<<" "<<featureArray.size()<<endl;
  ICP(featureArray, modelArray, /*settings,*/ currTranslation, currRotation, errorArray);

  cout<<"Translation="<<currTranslation<<endl;
  cout<<"Rotation="<<currRotation<<endl;

  return 0;
}


















