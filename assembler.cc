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

int ReadAssemblerConfigFile(string assemblerConfigFilename, struct AssemblerParams *assemblerParams)
{
  ifstream configFile (assemblerConfigFilename.c_str());
  std::string line;
  double val, val1, val2, val3, val4; 
  std::string identifier;

  if (configFile.is_open()){ 
 
    //matching params
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline0; 
    sline0 << line;
    sline0 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->matchingMode = val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline1; 
    sline1 << line;
    sline1 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->maxNumStarts= val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline2; 
    sline2 << line;
    sline2 >> identifier >> val1 >> val2;
    //cout<<val<<endl;
    assemblerParams->deltaLonLat(0)=val1;
    assemblerParams->deltaLonLat(1)=val2;
    
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline3; 
    sline3 << line;
    sline3 >> identifier >> val1 >> val2;
    //cout<<val<<endl;
    assemblerParams->samplingStep(0)=val1;
    assemblerParams->samplingStep(1)=val2;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline4; 
    sline4 << line;
    sline4 >> identifier >> val1 >> val2;
    //cout<<val<<endl;
    assemblerParams->matchWindowHalfSize(0)=val1;
    assemblerParams->matchWindowHalfSize(1)=val2;
    
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline5; 
    sline5 << line;
    sline5 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->maxNumIter;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline6; 
    sline6 << line;
    sline6 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->minConvThresh;
  
    //tiling params
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline10; 
    sline10<<line;
    sline10 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->tileSizeDEM=val;
    
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline11; 
    sline11<<line;
    sline11 >> identifier >> val1>>val2>>val3>>val4;
    //cout<<val<<endl;
    assemblerParams->paddingParamsDEM(0)=val1;
    assemblerParams->paddingParamsDEM(1)=val2;
    assemblerParams->paddingParamsDEM(2)=val3;
    assemblerParams->paddingParamsDEM(3)=val4;
    
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline12; 
    sline12<<line;
    sline12 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->tileSizeDRG=val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline13; 
    sline13<<line;
    sline13 >> identifier >> val1>>val2>>val3>>val4;
    //cout<<val<<endl;
    assemblerParams->paddingParamsDRG(0)=val1;
    assemblerParams->paddingParamsDRG(1)=val2;
    assemblerParams->paddingParamsDRG(2)=val3;
    assemblerParams->paddingParamsDRG(3)=val4;
    
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline14; 
    sline14<<line;
    sline14 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->foreNoDataValDEM = val;
  
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline15; 
    sline15<<line;
    sline15>> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->backNoDataValDEM=val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline16; 
    sline16<<line;
    sline16>> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->foreNoDataValDRG = val;
  
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline17; 
    sline17<<line;
    sline17>> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->backNoDataValDRG=val;
   
    configFile.close();
    return 1;
  }
  else{
    cout << "Unable to open settings file. Use default values"; 

     //matching params - START
     assemblerParams->matchingMode = 1; 
     assemblerParams->maxNumStarts = 36;
     assemblerParams->deltaLonLat(0) = 0.0001;
     assemblerParams->deltaLonLat(1) = 0.0001;
     assemblerParams->samplingStep(0) = 256;
     assemblerParams->samplingStep(1) = 256;
     assemblerParams->matchWindowHalfSize(0) = 9;
     assemblerParams->matchWindowHalfSize(1) = 9;
     assemblerParams->maxNumIter = 5;
     assemblerParams->minConvThresh = 0.01;
     //matching params - END

     //tiling params - START
     assemblerParams->tileSizeDEM=128;
     assemblerParams->paddingParamsDEM(0)=0;
     assemblerParams->paddingParamsDEM(1)=0;
     assemblerParams->paddingParamsDEM(2)=1;
     assemblerParams->paddingParamsDEM(3)=1;
     
     assemblerParams->tileSizeDRG=512;
     assemblerParams->paddingParamsDRG(0)=1;
     assemblerParams->paddingParamsDRG(1)=1;
     assemblerParams->paddingParamsDRG(2)=2;
     assemblerParams->paddingParamsDRG(3)=2;
     
     assemblerParams->foreNoDataValDEM = -3.4028226550889e+38;
     assemblerParams->backNoDataValDEM = -3.4028226550889e+38;
     assemblerParams->foreNoDataValDRG = 0;
     assemblerParams->backNoDataValDRG = 0;
     //matching params - END
  
     return 0;
  }
  
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
 
  struct AssemblerParams assemblerParams;
  string assemblerConfigFilename;
 
  if( ReadAssemblerConfigFile(assemblerConfigFilename, &assemblerParams) ){
     std::cerr << "Assembler Config file " << assemblerConfigFilename << " found." << endl;
  }
  else{
    std::cerr << "Assembler Config file " << assemblerConfigFilename << " not found, using defaults." << endl;
  }
  
  //PrintGlobalParams(&settings);
  //std::cerr << assemblerParams << endl;


  //determine the noDataValues for fore and back files - START
  float backNoDataVal;
  boost::shared_ptr<DiskImageResource> back_rsrc( new DiskImageResourceGDAL(backFile) );
  if (back_rsrc->has_nodata_read()){
      backNoDataVal = back_rsrc->nodata_read();
      cout<<"noDataVal for background ="<<backNoDataVal;
  }
  else{
      backNoDataVal = assemblerParams.backNoDataValDEM;
      cout<<"noDataVal not found for background, using: "<<backNoDataVal<<endl;
  }
    
  float foreNoDataVal;// = settings.noDataVal;
  boost::shared_ptr<DiskImageResource> fore_rsrc( new DiskImageResourceGDAL(foreFile) );
  if (fore_rsrc->has_nodata_read()){
    foreNoDataVal = fore_rsrc->nodata_read();
    cout<<"noDataVal for foreground ="<<foreNoDataVal;
  }
  else{
      foreNoDataVal = assemblerParams.foreNoDataValDEM;
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
   
    if (assemblerParams.matchingMode != 0){
       
       float matchError;
       Vector2 delta_lonlat;        
 
       int numVerRestarts = sqrt(assemblerParams.maxNumStarts); 
       int numHorRestarts = sqrt(assemblerParams.maxNumStarts);

       for (int k = -(numVerRestarts-1)/2; k < (numVerRestarts+1)/2; k++){
	 delta_lonlat(0) = k*assemblerParams.deltaLonLat(0); //~5 meters increments
	 for (int l = -(numHorRestarts-1)/2; l < (numHorRestarts+1)/2; l++){
	   delta_lonlat(1) = l*assemblerParams.deltaLonLat(1); //~5 meters increments
	   cout<<"k ="<<k<<", l="<<l<<endl; 
	   
	   cout<<"Feature extraction ..."<<endl;
	   vector<Vector3> featureArray = GetFeatures(foreDEM, foreDEMGeo, backDEM, backDEMGeo, 
						      assemblerParams.samplingStep, delta_lonlat, 
                                                      assemblerParams.foreNoDataValDEM);
	   cout<<"done."<<endl;
	   vector<float> errorArray;
	   errorArray.resize(featureArray.size());
	   
	   Vector3 currTranslation;
	   Matrix<float, 3,3 > currRotation;
	   Vector3 currCenter;
	   
	   cout<<"running ICP"<<endl;
           struct CoregistrationParams coregistrationParams;
           coregistrationParams.matchingMode = assemblerParams.matchingMode;
           coregistrationParams.reflectanceType = 0;
           coregistrationParams.analyseFlag = 0;
           coregistrationParams.useReflectanceFeatures = 1;
           coregistrationParams.topPercentFeatures = 1;
           coregistrationParams.samplingStep = assemblerParams.samplingStep;
           coregistrationParams.matchWindowHalfSize = assemblerParams.matchWindowHalfSize;
           coregistrationParams.maxNumIter = assemblerParams.maxNumIter;
           coregistrationParams.maxNumStarts = assemblerParams.maxNumStarts;
           coregistrationParams.noDataVal = assemblerParams.foreNoDataValDEM;
           coregistrationParams.minConvThresh = assemblerParams.minConvThresh;

	   ICP_DEM_2_DEM(featureArray, backDEM, backDEMGeo, foreDEMGeo, coregistrationParams,
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
  if ((mode.compare("DRG")==0) || (mode.compare("DEM_DRG")==0)){
  
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
 
}


















