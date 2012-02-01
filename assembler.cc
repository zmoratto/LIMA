// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>
#include <sstream>

//boost
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>

#include "coregister.h"
#include "icp.h"
#include "assembler.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

int main( int argc, char *argv[] ) {
   
  std::string foreFile;
  std::string backFile;
  std::string resDir;
  std::string configFilename = "assembler_settings.txt";
  std::string mode;//it can be DEM or DRG
  po::options_description general_options("Options");
  
  general_options.add_options()
    ("mode,m", po::value<std::string>(&mode)->default_value("DEM"), "mode.")
    ("backFile,b", po::value<std::string>(&backFile)->default_value(""), "background file.")
    ("resDir,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
    ("configFilename,c", po::value<std::string>(&configFilename)->default_value("assembler_settings.txt"), "configuration filename.")
    ("help,h", "Display this help message");
  
  po::options_description hidden_options("");
  
  hidden_options.add_options()
  ("foreFile, f", po::value<std::string> (&foreFile));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);
  
  po::positional_options_description p;
  p.add("foreFile", -1);
 
  std::ostringstream usage;
  usage << "Description: main code for DEM to DEM alignment" << std::endl << std::endl;
  usage << general_options << std::endl;
  
  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch ( po::error &e ) {
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }

  if( vm.count("help") ) {
    std::cerr << usage.str() << std::endl;
    return 1;
  }
    
  if( vm.count("foreFile") < 1 ) {
    std::cerr << "Error: Must specify at least one foreground file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }

  if (mode.compare("DEM")==0){
      printf("DEM\n");
  }
  else{
      printf("DRG\n");
  }
 
  
  struct CoregistrationParams settings;
  if( ReadConfigFile(configFilename, &settings) ){
     std::cerr << "Config file " << configFilename << " found." << endl;
  }
  else{
    std::cerr << "Config file " << configFilename << " not found, using defaults." << endl;
  }
  //PrintGlobalParams(&settings);
  std::cerr << settings << endl;
  
  Vector3 translation;
  Matrix<float, 3,3 > rotation;
  rotation[0][0]=1.0;
  rotation[1][1]=1.0;
  rotation[2][2]=1.0;

  Vector3 center;
  Vector2 bestDeltaLonLat;
  
  bestDeltaLonLat(0)=0;
  bestDeltaLonLat(1)=0;

  if ((mode.compare("DEM")==0) || (mode.compare("DEM_DRG")==0) ){
    string backDEMFilename = backFile;
    string foreDEMFilename = foreFile;
    //string assembledDEMFilename = resDir+"/assembled_dem.tif";
   
    //large image low res - background
    cout<<"opening"<<backDEMFilename<<"..."<<endl; 
    DiskImageView<float>   backDEM(backDEMFilename);
    GeoReference backDEMGeo;
    read_georeference(backDEMGeo, backDEMFilename);
        
    //small image high res - foreground
    cout<<"opening"<<foreDEMFilename<<"..."<<endl;  
    DiskImageView<float>  foreDEM(foreDEMFilename);
    cout<<"done."<<endl;

    cout<<"opening georeference for "<<foreDEMFilename<<"..."<<endl;  
    GeoReference foreDEMGeo;
    read_georeference(foreDEMGeo, foreDEMFilename);
    cout<<"done"<<endl;
 
    float minMatchError = 100000000.0;
    
    if (settings.matchingMode != 0){
       
       Vector2 delta_lonlat; 
       
       int numVerRestarts = sqrt(settings.maxNumStarts); 
       int numHorRestarts = sqrt(settings.maxNumStarts);
      
       for (int k = -(numVerRestarts-1)/2; k < (numVerRestarts+1)/2; k++){
	   delta_lonlat(0) = k*0.001; 
	   for (int l = -(numHorRestarts-1)/2; l < (numHorRestarts+1)/2; l++){
	       delta_lonlat(1) = l*0.001;
               cout<<"k ="<<k<<", l="<<l<<endl; 
	       cout<<"Feature extraction ..."<<endl;
	       vector<Vector3> featureArray = GetFeatures(foreDEM, foreDEMGeo, backDEM, backDEMGeo, 
                                                          settings.samplingStep, delta_lonlat, settings.noDataVal);
	       cout<<"done."<<endl;
               vector<float> errorArray;
	       errorArray.resize(featureArray.size());
               
               Vector3 currTranslation;
               Matrix<float, 3,3 > currRotation;
	       cout<<"running ICP"<<endl;
               ICP_DEM_2_DEM(featureArray, backDEM, backDEMGeo, foreDEMGeo, settings,
		             currTranslation, currRotation, center, errorArray);
               cout<<"done."<<endl;
                
	       float matchError = 0;
	       for (int m = 0; m < errorArray.size(); m++){
		 matchError = matchError+errorArray[m];
	       }
	       matchError = matchError/errorArray.size();
               cout<<"currMatchError="<<matchError<<", minMatchError="<<minMatchError<<endl;

	       if (matchError < minMatchError){
		  minMatchError = matchError;
                  bestDeltaLonLat = delta_lonlat;
		  rotation = currRotation;
		  translation = currTranslation;
	       }
	       
	 }
       }
       
       //bestDeltaLonLat(0) = 0.002;
       //bestDeltaLonLat(1) = 0.002;
    }
    
    cout<<"minMatchError "<<minMatchError<<endl;
    cout<<"final Rotation matrix "<<rotation<<endl;
    cout<<"final translation vector "<<translation<<endl;
    cout<<"bestDeltaLonLat="<<bestDeltaLonLat<<endl;


    vector<struct tilingParams> tileParamsArray;
    int tileSize = 128;
    int numTiles = 128;
    Vector2 DEMOffset;
    DEMOffset(0) = (tileSize*numTiles - backDEM.cols())/2;
    DEMOffset(1) = (tileSize*numTiles - backDEM.rows())/2;
    cout<<"DEMOffset="<<DEMOffset<<endl;

    ComputeBoundariesDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo, tileSize, DEMOffset, tileParamsArray);
    for (int i= 0; i <tileParamsArray.size(); i++){
        stringstream ss;
        ss<<tileParamsArray[i].horTileIndex<<"_"<<tileParamsArray[i].verTileIndex;
        //ss<<tileParamsArray[i].verTileIndex;

	string assembledDEMFilename =  resDir+"/assembled_"+ss.str()+"_dem.tif";
        cout<<assembledDEMFilename<<endl;
	ComputeAssembledDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo, assembledDEMFilename, 
			    translation, rotation, center, bestDeltaLonLat, tileParamsArray[i]);   
    }
    //ComputeAssembledDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo, assembledDEMFilename, 
    //                    translation, rotation, center, bestDeltaLonLat, tileParams);    
  }

  //DRG assembler
  if (mode.compare("DRG")==0){
  
      string backDRGFilename = backFile;//"../MSLData/Mars/MER_HIRISE/PSP_001777_1650_1m_o-crop-geo.tif";
      string foreDRGFilename = foreFile;//"../MSLData/Mars/MER_HIRISE/Photo-mod.tif";
      int tileSize = 512;
      int numTiles = 128;
      Vector2 DRGOffset;
      Vector4 paddingParams;

      //small image high res - foreground 
      //DiskImageView<PixelRGB<uint8> >  backDRG(backDRGFilename);
      DiskImageView<PixelGray<uint8> >  backDRG(backDRGFilename);
      GeoReference backDRGGeo;
      read_georeference(backDRGGeo, backDRGFilename);
      printf("done opening the the backDRG\n");
      
      //large image low res - background
      DiskImageView<PixelRGB<uint8> >  foreDRG(foreDRGFilename);
      GeoReference foreDRGGeo;
      read_georeference(foreDRGGeo, foreDRGFilename);
      printf("done opening the the foreDRG\n");

      DRGOffset(0) = (tileSize*numTiles - backDRG.cols())/2;
      DRGOffset(1) = (tileSize*numTiles - backDRG.rows())/2;
      cout<<"DRGOffset="<<DRGOffset<<endl;
   
      vector<struct tilingParams> tileParamsArray;
      
      ComputeBoundariesDRG(foreDRG, foreDRGGeo, backDRG, backDRGGeo, tileSize, DRGOffset, tileParamsArray);
      for (int i= 0; i <tileParamsArray.size(); i++){
        stringstream ss;
        //ss<<i;
        ss<<tileParamsArray[i].horTileIndex<<"_"<<tileParamsArray[i].verTileIndex;
	string assembledDRGFilename =  resDir+"/assembled_"+ss.str()+"_drg.tif";
        cout<<assembledDRGFilename<<endl;
	ComputeAssembledDRG(foreDRG, foreDRGGeo, backDRG, backDRGGeo, assembledDRGFilename, 
			    translation, rotation, center, bestDeltaLonLat, tileParamsArray[i]);   
      }
      
   }

   if (mode.compare("DEM_DRG")==0){
      string backDRGFilename = "../../../msl/MSLData/Mars/MER_HIRISE/PSP_001777_1650_1m_o-crop-geo.tif";
      string foreDRGFilename = "../../../msl/MSLData/Mars/MER_HIRISE/Photo-mod.tif";
      string assembledDRGFilename =  resDir+"/assembled_drg.tif";
 
      //small image high res - foreground 
      DiskImageView<PixelGray<uint8> >  backDRG(backDRGFilename);
      GeoReference backDRGGeo;
      read_georeference(backDRGGeo, backDRGFilename);
      printf("done opening the the backDRG\n");

      //large image low res - background
      DiskImageView<PixelRGB<uint8> >  foreDRG(foreDRGFilename);
      GeoReference foreDRGGeo;
      read_georeference(foreDRGGeo, foreDRGFilename);
      printf("done opening the the foreDRG\n");
   }
    
}


















