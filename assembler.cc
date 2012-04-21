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


void SaveAssembledPC(string DEMFilename, string assembledPCFilename)
{
  cout<<"Saving the PC to file ..."<<endl;
  
  DiskImageView<float>   DEM(DEMFilename);
  GeoReference Geo;
  read_georeference(Geo, DEMFilename);
  int verStep = 20;
  int horStep = 20;
  float noDataVal = -37687.0;
  
 
  ImageViewRef<float>interpDEM = interpolate(edge_extend(DEM.impl(),
							 ConstantEdgeExtension()),
					     BilinearInterpolation());

  Vector2 pix_0(0,0);
  Vector2 lonlat_0= Geo.pixel_to_lonlat(pix_0);
  Vector3 lonlat3_0(lonlat_0(0), lonlat_0(1), (interpDEM.impl())(0,0));
  Vector3 xyz_0 = Geo.datum().geodetic_to_cartesian(lonlat3_0);
  
  FILE* fp = fopen(assembledPCFilename.c_str(), "wt");
   
  for (int j = 0; j < DEM.impl().rows(); j=j+verStep){
    for (int i = 0; i < DEM.impl().cols(); i=i+horStep){
      if ((isnan(DEM.impl()(i,j))!=FP_NAN) && (DEM.impl()(i,j)) != -noDataVal)  {
	Vector2 pix(i,j);
	Vector2 lonlat = Geo.pixel_to_lonlat(pix);
	Vector3 lonlat3(lonlat(0), lonlat(1), (interpDEM.impl())(i,j));
	Vector3 xyz = Geo.datum().geodetic_to_cartesian(lonlat3);
	fprintf(fp, "%f %f %f\n", xyz[0]-xyz_0[0], xyz[1]-xyz_0[1], xyz[2]-xyz_0[2]);
      }
    }
  }

  fclose(fp);
}
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
  
  struct AssemblerParams assemblerParams;
  assemblerParams.deltaLonLat(0)=0.0001;
  assemblerParams.deltaLonLat(1)=0.0001;
 
  assemblerParams.tileSizeDEM=128;
  assemblerParams.paddingParamsDEM(0)=0;
  assemblerParams.paddingParamsDEM(1)=0;
  assemblerParams.paddingParamsDEM(2)=1;
  assemblerParams.paddingParamsDEM(3)=1;

  assemblerParams.tileSizeDRG=512;
  assemblerParams.paddingParamsDRG(0)=1;
  assemblerParams.paddingParamsDRG(1)=1;
  assemblerParams.paddingParamsDRG(2)=2;
  assemblerParams.paddingParamsDRG(3)=2;

  assemblerParams.foreNoDataValDEM = -3.4028226550889e+38;
  assemblerParams.backNoDataValDEM = -3.4028226550889e+38;
  assemblerParams.foreNoDataValDRG = 0;
  assemblerParams.backNoDataValDRG = 0;
 
  assemblerParams.matchingMode = 0; 

  //determine the noDataValues for fore and back files - START
  float backNoDataVal = settings.noDataVal;
  boost::shared_ptr<DiskImageResource> back_rsrc( new DiskImageResourceGDAL(backFile) );
  if (back_rsrc->has_nodata_read()){
      backNoDataVal = back_rsrc->nodata_read();
      cout<<"noDataVal for background ="<<backNoDataVal;
  }
  else{
      backNoDataVal = assemblerParams.backNoDataValDEM;//settings.noDataVal;
      cout<<"noDataVal not found for background, using: "<<backNoDataVal<<endl;
  }
    
  float foreNoDataVal = settings.noDataVal;
  boost::shared_ptr<DiskImageResource> fore_rsrc( new DiskImageResourceGDAL(foreFile) );
  if (fore_rsrc->has_nodata_read()){
    foreNoDataVal = fore_rsrc->nodata_read();
    cout<<"noDataVal for foreground ="<<foreNoDataVal;
  }
  else{
      foreNoDataVal = assemblerParams.foreNoDataValDEM;//settings.noDataVal;
      cout<<"noDataVal not found for foreground, using: "<<foreNoDataVal<<endl;
  }
  //determine the noDataValues for fore and back files - END

  RegistrationParams registrationParams;
  
  registrationParams.translation(0)=0.0;
  registrationParams.translation(1)=0.0;
  registrationParams.translation(2)=0.0;
  
  registrationParams.rotation(0,0)=1.0;
  registrationParams.rotation(0,1)=0.0;
  registrationParams.rotation(0,2)=0.0;
  registrationParams.rotation(1,0)=0.0; 
  registrationParams.rotation(1,1)=1.0;   
  registrationParams.rotation(1,2)=0.0;
  registrationParams.rotation(2,0)=0.0;   
  registrationParams.rotation(2,1)=0.0;
  registrationParams.rotation(2,2)=1.0;

  registrationParams.center(0)=0.0;
  registrationParams.center(1)=0.0;
  registrationParams.center(2)=0.0;

  registrationParams.bestDeltaLonLat(0)=0.0;
  registrationParams.bestDeltaLonLat(1)=0.0;
  
  registrationParams.error=100000000.0;
  
  if ((mode.compare("DEM")==0) || (mode.compare("DEM_DRG")==0) ){

    string backDEMFilename = backFile;
    string foreDEMFilename = foreFile;
 
    //large image low res - background
    cout<<"opening"<<backDEMFilename<<"..."<<endl; 
    DiskImageView<float>   backDEM(backDEMFilename);
    GeoReference backDEMGeo;
    read_georeference(backDEMGeo, backDEMFilename);
    
    //small image high res - foreground
    cout<<"opening"<<foreDEMFilename<<"..."<<endl;  
    DiskImageView<float>  foreDEM(foreDEMFilename);
    GeoReference foreDEMGeo;
    read_georeference(foreDEMGeo, foreDEMFilename);
   
    if (settings.matchingMode != 0){
       
       float matchError;
       Vector2 delta_lonlat;        
       int numVerRestarts = sqrt(settings.maxNumStarts); 
       int numHorRestarts = sqrt(settings.maxNumStarts);
      
       for (int k = -(numVerRestarts-1)/2; k < (numVerRestarts+1)/2; k++){
	 delta_lonlat(0) = k*assemblerParams.deltaLonLat(0); //~5 meters increments
	 for (int l = -(numHorRestarts-1)/2; l < (numHorRestarts+1)/2; l++){
	   delta_lonlat(1) = l*assemblerParams.deltaLonLat(1); //~5 meters increments
	   cout<<"k ="<<k<<", l="<<l<<endl; 
	   
	   cout<<"Feature extraction ..."<<endl;
	   vector<Vector3> featureArray = GetFeatures(foreDEM, foreDEMGeo, backDEM, backDEMGeo, 
						      settings.samplingStep, delta_lonlat, settings.noDataVal);
	   cout<<"done."<<endl;
	   vector<float> errorArray;
	   errorArray.resize(featureArray.size());
	   
	   Vector3 currTranslation;
	   Matrix<float, 3,3 > currRotation;
	   Vector3 currCenter;
	   
	   cout<<"running ICP"<<endl;
	   ICP_DEM_2_DEM(featureArray, backDEM, backDEMGeo, foreDEMGeo, settings,
			 currTranslation, currRotation, currCenter, matchError);
	   cout<<"done."<<endl;
	   
	   if (matchError <registrationParams.error){       
	     registrationParams.error = matchError;
	     registrationParams.bestDeltaLonLat = delta_lonlat;
	     registrationParams.rotation = currRotation;
	     registrationParams.translation = currTranslation;
	     registrationParams.center = currCenter;
	   }   
	 }
       }
    }
  
    cout<<"minMatchError "<<registrationParams.error<<endl;
    cout<<"final Rotation matrix "<<registrationParams.rotation<<endl;
    cout<<"final translation vector "<<registrationParams.translation<<endl;
    cout<<"center="<<registrationParams.center<<endl;
    cout<<"bestDeltaLonLat="<<registrationParams.bestDeltaLonLat<<endl;

    vector<struct TilingParams> tileParamsArray;

    Vector4 lonlatBB;
    ComputeLonLatBoxDEM(foreDEM, foreDEMGeo, assemblerParams.foreNoDataValDEM, 
                        registrationParams.bestDeltaLonLat, lonlatBB);
    
    ComputeBoundaries(foreDEM, foreDEMGeo, backDEM, backDEMGeo,
                      registrationParams, assemblerParams, 
                      0, lonlatBB, tileParamsArray);
    
    for (int i= 0; i <tileParamsArray.size(); i++){
	ComputeAssembledDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo,
                            assemblerParams.foreNoDataValDEM, resDir,  
                            registrationParams, tileParamsArray[i]);   
    }
   
  }

  //DRG assembler
  if (mode.compare("DRG")==0){
  
      Vector2 DRGOffset;
      //small image high res - foreground 
      //DiskImageView<PixelRGB<uint8> >  backDRG(backFile);
      DiskImageView<PixelGray<uint8> >  backDRG(backFile);
      GeoReference backDRGGeo;
      read_georeference(backDRGGeo, backFile);
      printf("done opening the the backDRG\n");
      
      //large image low res - background
      DiskImageView<PixelRGB<uint8> >  foreDRG(foreFile);
      GeoReference foreDRGGeo;
      read_georeference(foreDRGGeo, foreFile);
      printf("done opening the the foreDRG\n");
   
      vector<struct TilingParams> tileParamsArray;
     
      //determines the background DEM tiles to be modified by the new position of the foreground DRG 
      Vector4 lonlatBB;
      ComputeLonLatBoxDRG(foreDRG, foreDRGGeo, assemblerParams.foreNoDataValDRG, 
                          registrationParams.bestDeltaLonLat, lonlatBB);      
      ComputeBoundaries(foreDRG, foreDRGGeo, backDRG, backDRGGeo, 
                        registrationParams, assemblerParams, 
                        1, lonlatBB, tileParamsArray);

      for (int i= 0; i <tileParamsArray.size(); i++){
     
	ComputeAssembledDRG(foreDRG, foreDRGGeo, backDRG, backDRGGeo, 
                            assemblerParams.foreNoDataValDRG, resDir, 
                            registrationParams, tileParamsArray[i]);   
      }
      
   }

   // this will most likely no longer be used - START
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
    // this will most likely no longer be used - END  
  
}


















