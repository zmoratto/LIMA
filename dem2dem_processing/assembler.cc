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

#include "icp.h"
#include "assembler.h"
#include "pyr_tiling.h"
#include "../common/string_util.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;


float ComputeDEMAccuracy(GeoReference foreGeo, Vector2 forePix, float backAccuracy)
{
   
  Vector2 forePoint;
  forePoint = foreGeo.pixel_to_point(forePix);
  
  float distToCam = sqrt(forePoint(0)*forePoint(0) + forePoint(1)*forePoint(1));
  float foreWeight = 0.5;
  
  //delta_d = 43 micro, f = 43mm, b = 0.3m
  float foreAccuracy = distToCam*distToCam/300.0; ///delta_r = r*r*delta_d/(b*f)
  foreWeight = backAccuracy/(foreAccuracy+backAccuracy);
  
  return foreWeight;

}
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
   
    cout << "Reading settings file "<< assemblerConfigFilename<<"."<<endl; 
   
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
    stringstream sline01; 
    sline01 << line;
    sline01 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->weightingMode = val;

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
    stringstream sline200; 
    sline200 << line;
    sline200 >> identifier >> val1;
    //cout<<val<<endl;
    assemblerParams->useForeLonLatRadOffset=val1;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline20; 
    sline20 << line;
    sline20 >> identifier >> val1 >> val2;
    //cout<<val<<endl;
    assemblerParams->foreLonLatOrigin(0)=val1;
    assemblerParams->foreLonLatOrigin(1)=val2;
    //assemblerParams->foreLonLatOrigin(2)=val3;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline201; 
    sline201 << line;
    sline201 >> identifier >> val1;
    //cout<<val<<endl;
    assemblerParams->foreMaxPPD=val1;


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
    assemblerParams->maxNumIter=val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline6; 
    sline6 << line;
    sline6 >> identifier >> val;
    //cout<<val<<endl;
    assemblerParams->minConvThresh=val;
  
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
    cout << "Unable to open settings file "<< assemblerConfigFilename<<". Use default values"<<endl; 

     //matching params - START
     assemblerParams->matchingMode = 1; 
     assemblerParams->weightingMode = 0;
     assemblerParams->maxNumStarts = 36;
     assemblerParams->deltaLonLat(0) = 0.0001;
     assemblerParams->deltaLonLat(1) = 0.0001;
     assemblerParams->foreLonLatOrigin(0) = 0.0;//137.4-180
     assemblerParams->foreLonLatOrigin(1) = 0.0;//-4.5
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
  std::string backDRGFile;
  std::string resDir;
  std::string configFilename = "assembler_settings.txt";
  std::string mode;//it can be DEM or DRG
  po::options_description general_options("Options");
  
  general_options.add_options()
    ("mode,m", po::value<std::string>(&mode)->default_value("DEM"), "mode.")
    ("backFile,b", po::value<std::string>(&backFile)->default_value(""), "background file.")
    ("backImgFile,p", po::value<std::string>(&backDRGFile)->default_value(""), "background photo file.")
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
  ReadAssemblerConfigFile(configFilename, &assemblerParams);

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
  registrationParams.deltaRad=0.0;
  
  registrationParams.error=std::numeric_limits<float>::max();

  cout<<"foreDEMFilename="<<foreFile<<endl;
  cout<<"backDEMFilename="<<backFile<<endl;
  cout<<"backDRGFilename="<<backDRGFile<<endl;

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

    /*
    cout.precision(9);
    Vector2 lonlat_temp;
    lonlat_temp(0) = 137.4417;
    lonlat_temp(1) = -4.5895;
    Vector2 point = backDEMGeo.lonlat_to_point(lonlat_temp);
    cout<<"POINT="<<point<<endl;
    Vector2 pixel = backDEMGeo.lonlat_to_pixel(lonlat_temp);
    cout<<"PIXEL="<<pixel<<endl;
    cout<<"height="<<backDEM(pixel(0), pixel(1));
    */

    vector<Vector3> featureArray;

    Vector2 camPoint;
    camPoint(0) = 0;
    camPoint(1) = 0;
    Vector2 camLonLat = foreDEMGeo.point_to_lonlat(camPoint);
    cout<<"camLonLat="<<camLonLat<<endl;

    
    //assemblerParams.foreLonLatOrigin = backDEMGeo.lonlat_to_point(assemblerParams.foreLonLatOrigin);
    cout<<"fore_point_orig="<< assemblerParams.foreLonLatOrigin<<endl;

    //pixel to lonlat of the initial point - START
    assemblerParams.foreLonLatOrigin = backDEMGeo.point_to_lonlat(assemblerParams.foreLonLatOrigin);
    cout<<"fore_lonlat_orig="<< assemblerParams.foreLonLatOrigin<<endl;
    //pixel to lonlat of the initial point - END

 
    Vector2 offsetLonLat;
    //compute the lonlat offset between the original lonlat origin of the foregoround and the new origin in the configuration file.  
    if (assemblerParams.useForeLonLatRadOffset == 0){
      //use current values of the foreground origin
      offsetLonLat(0) = 0.0;
      offsetLonLat(1) = 0.0;
    }
    else{//use new values of the foregoround origin
      offsetLonLat(0) = assemblerParams.foreLonLatOrigin(0) - camLonLat(0);    
      offsetLonLat(1) = assemblerParams.foreLonLatOrigin(1) - camLonLat(1);
    }
    cout<<"offsetLonLat="<<offsetLonLat<<endl; 
    
    //matchingMode = 0 do nothing, no alignmment
    //matchingMode = 1 altitude alignment, no ICP
    //matchingMode = 2 ICP, translation only
    //matchingMode = 3 ICP, rotation and translation
 
    if (assemblerParams.matchingMode == 1){//altitude alignment
      cout<<"altitude only alignment"<<endl;
      float deltaRad = ComputeAverageDiff(foreDEM, foreDEMGeo, backDEM, backDEMGeo, assemblerParams.foreNoDataValDEM, offsetLonLat, registrationParams);
      registrationParams.deltaRad = deltaRad;
      registrationParams.bestDeltaLonLat = offsetLonLat;
    }

    if (assemblerParams.matchingMode > 1){//run ICP

       cout<<"ICP alignment"<<endl;

       float matchError;
       Vector2 delta_lonlat;        
 
       int numVerRestarts = sqrt(assemblerParams.maxNumStarts); 
       int numHorRestarts = sqrt(assemblerParams.maxNumStarts);  
      
       for (int k = -(numVerRestarts-1)/2; k < (numVerRestarts+1)/2; k++){
	// delta_lonlat(0) = offsetLonLat(0) +k*assemblerParams.deltaLonLat(0); //~5 meters increments
        delta_lonlat(0) = k*assemblerParams.deltaLonLat(0); //~5 meters increments
	 for (int l = -(numHorRestarts-1)/2; l < (numHorRestarts+1)/2; l++){
	   //delta_lonlat(1) = offsetLonLat(1) +assemblerParams.deltaLonLat(1); //~5 meters increments
	   delta_lonlat(1) = l*assemblerParams.deltaLonLat(1); //~5 meters increments
	   cout<<"k ="<<k<<", l="<<l<<", delta_lon="<<delta_lonlat(0)<<", delta_lat="<<delta_lonlat(1)<<endl; 
	   cout<<"Foreground DEM feature extraction ..."<<endl;
	
	   featureArray = GetFeatures(foreDEM, foreDEMGeo,  
				      assemblerParams.samplingStep, offsetLonLat + delta_lonlat, 
				      assemblerParams.foreNoDataValDEM);
	   
	   cout<<"done."<<endl;
	   vector<float> errorArray;
	   errorArray.resize(featureArray.size());
	   
	   Vector3 currTranslation = registrationParams.translation;
	   Matrix<float, 3,3 > currRotation = registrationParams.rotation;
	   Vector3 currCenter = registrationParams.center;
	   
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
           cout<<"coregParams="<<coregistrationParams<<endl;

         
	   ICP_DEM_2_DEM(featureArray, backDEM, backDEMGeo, foreDEMGeo, assemblerParams.matchingMode, 
			 coregistrationParams, currTranslation, currRotation, currCenter, matchError);
	   
	   cout<<"currCenter:"<<currCenter<<endl;
	   cout<<"currTranslation:"<<currTranslation<<endl;
           cout<<"currRotation:"<<currRotation<<endl;
           cout<<"------------------"<<endl;

	   /*
	   // make  a copy the featureArray for debug purposes
	   vector<Vector3> transfFeatureArray;
           transfFeatureArray.resize(featureArray.size());
	  
           for (int m = 0; m < featureArray.size(); m++){
	     transfFeatureArray[m](0) = featureArray[m](0);
             transfFeatureArray[m](1) = featureArray[m](1);
	     transfFeatureArray[m](2) = featureArray[m](2);
           }

           TransformFeatures( transfFeatureArray, currCenter, currTranslation, currRotation );     

        
	   
	   for (int m = 0; m < featureArray.size(); m++){
	     Vector3 transf_fore_lon_lat_rad = foreDEMGeo.datum().cartesian_to_geodetic(transfFeatureArray[m]);
	     Vector3 origin_fore_lon_lat_rad = foreDEMGeo.datum().cartesian_to_geodetic(featureArray[m]);
	     cout<<"transf_feature :"<<transfFeatureArray[m]<<", orig:"<<featureArray[m]<<endl;
             cout<<"transf_feature :"<<transf_fore_lon_lat_rad<<", orig:"<<origin_fore_lon_lat_rad<<endl;
	   }
           transfFeatureArray.clear();
	   */

           if (matchError <registrationParams.error){     
             cout<<"using initial registration parameters."<<endl;
	     registrationParams.error = matchError;
	     registrationParams.bestDeltaLonLat = offsetLonLat + delta_lonlat;
	     registrationParams.rotation = currRotation;
	     registrationParams.translation = currTranslation;
	     registrationParams.center = currCenter;
	   }
   
	 }//l
       }//k

       cout<<"done with ICP"<<endl;
       if (registrationParams.error ==std::numeric_limits<float>::max()){
	 cout<<"The background and foreground terrains cannot be matched"<<endl;
         
	 return 0;
       }

       cout<<"minMatchError "<<registrationParams.error<<endl;
       cout<<"final Rotation matrix "<<registrationParams.rotation<<endl;
       cout<<"final translation vector "<<registrationParams.translation<<endl;
       cout<<"center="<<registrationParams.center<<endl;
       cout<<"bestDeltaLonLat="<<registrationParams.bestDeltaLonLat<<endl;
       cout<<"cleanBestLonLat = "<<registrationParams.bestDeltaLonLat-offsetLonLat<<endl;

    }//done with ICP mode
  
   
    //TODO: need to create here a unified assembled georeference
    //it should be registered to the DEM
    //it will be passed to all the functions below instead of backDEMGeo.
    //this will attempt to remove the offset between DEM and DRG tiles.

    vector<struct TilingParams> tileParamsArray;
    Vector4 foreLonLatBB;

    cout<<"ASSEMBLER: ComputeLonLatBoxDEM of the foreground region"<<endl;
    ComputeLonLatBoxDEM(foreDEM, foreDEMGeo, assemblerParams.foreNoDataValDEM, 
                        registrationParams.bestDeltaLonLat, registrationParams, 
			foreLonLatBB);

    cout<<"ASSEMBLER: DEM bounding box: "<<foreLonLatBB<<endl;
    cout<<"ASSEMBLER: ComputeTileParams for DEMs"<<endl;

    int backDEMWidth = backDEM.cols(); 
    int backDEMHeight = backDEM.rows();
    
    Matrix<double> H;
    H = backDEMGeo.transform();
    cout<<"DEMGeo matrix="<<H<<endl;

    ComputeTileParams(backDEMWidth, backDEMHeight, foreDEMGeo, backDEMGeo,
                      assemblerParams.tileSizeDEM, assemblerParams.paddingParamsDEM, 
                      assemblerParams.foreNoDataValDEM, 0, 
                      foreLonLatBB, tileParamsArray);
   
    cout<<"ASSEMBLER: COMPUTE_ASSEMBLED_DEM for each tile"<<endl;
    
    for (int i= 0; i <tileParamsArray.size(); i++){
      
      ComputeAssembledDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo,
			  assemblerParams.foreNoDataValDEM, assemblerParams.foreLonLatOrigin,
                          assemblerParams.weightingMode, 
			  resDir, registrationParams, tileParamsArray[i]);   
      
    }
    
    if (mode.compare("DEM_DRG")==0){
  
      //large image low res - background
      cout<<"ASSEMBLER: backDRG filename="<<backDRGFile<<endl;

      //DiskImageView<PixelRGB<uint8> >  backDRG(backFile);
      DiskImageView<PixelGray<uint8> >  backDRG(backDRGFile);
      GeoReference backDRGGeo;  
      read_georeference(backDRGGeo, backDRGFile);
      
      //small image high res - foreground 
      string foreDRGFile = foreDEMFilename;
      cout<<"ASSEMBLER: foreDRG filename="<<foreDRGFile<<endl;
      FindAndReplace( foreDRGFile, string("Height"), string("Photo"));
   
      DiskImageView<PixelRGB<uint8> >  foreDRG(foreDRGFile);
      GeoReference foreDRGGeo;
      read_georeference(foreDRGGeo, foreDRGFile);

         
      //set the assembledDRGGeo NEW AND RISKY  - START
      //creates a unified assembled georeference for DRG and DEM
      //creates an assembledDRGGeo that is a copy of the backDEMGeo (or assembledDEMGeo)
      //except that the pixel scale is 4 times smaller in each dimension.
      //it creates the image size that is 4x4 bigger than the backDEM (or assembledDEM) image.
     
      
      Matrix<double> H;
      H = backDEMGeo.transform();
      cout<<"DEMGeo matrix="<<H<<endl;
     
      //H(0,2) = H(0,2)+(3.0/8.0)*H(0,0);
      //H(1,2) = H(1,2)+(3.0/8.0)*H(1,1);
      H(0,2) = H(0,2)+0.5*H(0,0);
      H(1,2) = H(1,2)+0.5*H(1,1);
      H(0,0) /=4.0;
      H(1,1) /=4.0;

      GeoReference assembledDRGGeo = backDEMGeo; 
      assembledDRGGeo.set_transform(H);
      cout<<"assembledDRGGeo matrix="<<H<<endl;
      H = backDEMGeo.transform();
      cout<<"DEMGeo matrix="<<H<<endl;
      
      /*
      GeoReference assembledDRGGeo = resample(backDEMGeo, 4.0, 4.0);
      Matrix<double> H;
      H = assembledDRGGeo.transform();
      cout<<"assembledDRGGeo matrix="<<H<<endl;
      */
      int assembledDRGWidth = backDEM.cols()*4;
      int assembledDRGHeight = backDEM.rows()*4;
      //set the assembledDRGGeo  NEW AND RISKY - END

      cout<<"ASSEMBLER: DRG bounding box: "<<foreLonLatBB<<endl;
      cout<<"ASSEMBLER: ComputeTileParams for DRG"<<endl;

      ComputeTileParams(assembledDRGWidth, assembledDRGHeight, foreDRGGeo, /*backDRGGeo*/assembledDRGGeo, 
			assemblerParams.tileSizeDRG, assemblerParams.paddingParamsDRG, 
                        assemblerParams.foreNoDataValDRG, 1, 
                        foreLonLatBB, tileParamsArray);


  
      /*
      //old style works well except the misalignment of DEM and DRG
      int backDRGWidth = backDRG.cols(); 
      int backDRGHeight = backDRG.rows();
      ComputeTileParams(backDRGWidth, backDRGHeight, foreDRGGeo, backDRGGeo, 
			assemblerParams.tileSizeDRG, assemblerParams.paddingParamsDRG, 
                        assemblerParams.foreNoDataValDRG, 1, 
                        foreLonLatBB, tileParamsArray);
      */
      cout<<"ASSEMBLER: ComputeAssembledDRG for each tile"<<endl;
     

      //assembledDRGGeo is the original backDEMGeo*4
      for (int i= 0; i <tileParamsArray.size(); i++){   
	
	//correct interpolation of the data, correct position.
        ComputeAssembledDRG(foreDRG, foreDRGGeo, 
                            backDRG, backDRGGeo,
                            backDEM, backDEMGeo,
                            assembledDRGGeo,
                            assemblerParams.foreNoDataValDRG, resDir, 
                            registrationParams, tileParamsArray[i]);
      }
    }
  }

  //DRG-only assembler does not support alignment options.
  //not that interesting for now
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
   
      Vector2 camPoint;
      camPoint(0) = 0;
      camPoint(1) = 0;
      Vector2 camLonLat = foreDRGGeo.point_to_lonlat(camPoint);
      cout<<"camLonLat="<<camLonLat<<endl;
      Vector2 offsetLonLat;
      if ((assemblerParams.foreLonLatOrigin(0) == 0) && (assemblerParams.foreLonLatOrigin(1) == 0)){
	offsetLonLat(0) = 0.0;
	offsetLonLat(1) = 0.0;
      }
      else{
	offsetLonLat(0) = assemblerParams.foreLonLatOrigin(0) - camLonLat(0);    
	offsetLonLat(1) = assemblerParams.foreLonLatOrigin(1) - camLonLat(1);
      }
      cout<<"offsetLonLat="<<offsetLonLat<<endl;
      

      vector<struct TilingParams> tileParamsArray;
      
      //determines the background DEM tiles to be modified by the new position of the foreground DRG 
      Vector4 foreLonLatBB;
      cout<<"computing foreLonLatBB...";

      registrationParams.bestDeltaLonLat = registrationParams.bestDeltaLonLat + offsetLonLat;

      cout<<"bestDeltaLonLat="<<registrationParams.bestDeltaLonLat<<endl;

      ComputeLonLatBoxDRG(foreDRG, foreDRGGeo, assemblerParams.foreNoDataValDRG, 
                          registrationParams.bestDeltaLonLat, foreLonLatBB);   
      
      cout<<"foreLonLatBB="<<foreLonLatBB<<endl;
      //cout<<"must be: 137.4,137.4,-4.50007,-4.49985"<<endl;
      //foreLonLatBB=Vector4(137.4,137.4,-4.50007,-4.49985)

      cout<<"computing boundaries...";

      int backDRGWidth = backDRG.cols(); 
      int backDRGHeight = backDRG.rows();
      ComputeTileParams(backDRGWidth, backDRGHeight, foreDRGGeo, backDRGGeo, 
			//assemblerParams,
			assemblerParams.tileSizeDRG, assemblerParams.paddingParamsDRG, assemblerParams.foreNoDataValDRG,   
			1, foreLonLatBB, tileParamsArray);
      cout<<"done";
      
      
      for (int i= 0; i <tileParamsArray.size(); i++){
	/*
        //this function does not alignment of the two images.
	ComputeAssembledDRG_old(foreDRG, foreDRGGeo, backDRG, backDRGGeo, 
                                assemblerParams.foreNoDataValDRG, resDir, 
                                registrationParams, tileParamsArray[i]);   
	*/
	
        ComputeAssembledDRG(foreDRG, foreDRGGeo, 
                            backDRG, backDRGGeo,
                            assemblerParams.foreNoDataValDRG, resDir,
                            registrationParams, tileParamsArray[i]);
	
      }
      
   }
 
}


















