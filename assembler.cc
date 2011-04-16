// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
#include <vw/Math.h>

#include "coregister.h"
#include "io.h"
#include "icp.h"
#include "assembler.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;

float ComputePixPerDegree(GeoReference Geo, int width, int height, int useUSGS_lonlat)
{

  //cout << Geo <<"\n";
  float radius = Geo.datum().semi_major_axis();
  cout<<"radius="<<radius<<endl;

  printf("width = %d, height = %d\n", width, height);
  Vector2 leftTopPixel(0,0);
  Vector2 leftTopLonLat = Geo.pixel_to_lonlat(leftTopPixel);


  Vector2 rightBottomPixel(width-1, height-1);
  Vector2 rightBottomLonLat = Geo.pixel_to_lonlat(rightBottomPixel);
  
  float minLon = leftTopLonLat(0);
  float minLat = leftTopLonLat(1);
  float maxLon = rightBottomLonLat(0);
  float maxLat = rightBottomLonLat(1);

  if (useUSGS_lonlat == 1){
    float usgs_2_lonlat = 180/(3.14159265*3396190);
     minLon = minLon*usgs_2_lonlat; 
     maxLon = maxLon*usgs_2_lonlat; 
     minLat = minLat*usgs_2_lonlat; 
     maxLat = maxLat*usgs_2_lonlat;
  }

  if (minLat < 0) minLat = minLat+180;
  if (minLon < 0) minLon = minLon+180;
  if (maxLat < 0) maxLat = maxLat+180;
  if (maxLon < 0) maxLon = maxLon+180;

  printf("minLon=%f, maxLon=%f, minLat=%f, maxLat=%f\n",
  	  minLon, maxLon, minLat, maxLat);

  float numPixPerDegree;
  numPixPerDegree = width/fabs(maxLon-minLon);
  numPixPerDegree = height/fabs(maxLat-minLat);
  
  return numPixPerDegree;
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
  usage << "Description: main code for DRG and DEM alignment" << std::endl << std::endl;
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
    std::cerr << "Error: Must specify at least one orthoprojected image file!" << std::endl << std::endl;
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
  ReadConfigFile((char*)configFilename.c_str(), &settings);
  PrintGlobalParams(&settings);
  
  Vector3 translation;
  Matrix<float, 3,3 > rotation;
  rotation[0][0]=1.0;
  rotation[1][1]=1.0;
  rotation[2][2]=1.0;

  if (mode.compare("DEM")==0){
    string backDEMFilename = backFile;
    string foreDEMFilename = foreFile;
    string assembledDEMFilename = resDir+"/assembled_dem.tif";
  
    //large image low res - background
    cout<<"opening"<<backDEMFilename<<"..."<<endl; 
    DiskImageView<PixelGray<float> >  backDEM(backDEMFilename);
    GeoReference backDEMGeo;
    read_georeference(backDEMGeo, backDEMFilename);
    float back_radius = backDEMGeo.datum().semi_major_axis();
    cout<<"radius="<<back_radius<<endl;
    cout<<"done"<<endl;
   
    //small image high res - foreground
    cout<<"opening"<<foreDEMFilename<<"..."<<endl;  
    DiskImageView<PixelGray<float> >  foreDEM(foreDEMFilename);
    GeoReference foreDEMGeo;
    read_georeference(foreDEMGeo, foreDEMFilename);
    float fore_radius = foreDEMGeo.datum().semi_major_axis();
    cout<<"radius="<<fore_radius<<endl;
    cout<<"done"<<endl;
   
    Vector2 bestDeltaLonLat;
    bestDeltaLonLat(0)=0;
    bestDeltaLonLat(1)=0;
    float minMatchError = 100000000.0;

    if (settings.matchingMode != 0){
       
       Vector2 delta_lonlat; 
      
       for (int k = -2; k < 3; k++){
	   delta_lonlat(0) = k*0.001; 
	   for (int l = -2; l < 3; l++){
	       delta_lonlat(1) = l*0.001;
	       printf("feature extraction ...\n");
	    
	       vector<Vector3> featureArray = GetFeatures(foreDEM, foreDEMGeo, backDEM, backDEMGeo, 
                                                          settings.samplingStep, delta_lonlat);
	       vector<float> errorArray;
	       errorArray.resize(featureArray.size());
               
               Vector3 currTranslation;
               Matrix<float, 3,3 > currRotation;
	       RunICP(featureArray, backDEM, backDEMGeo, foreDEMGeo, settings,
		      currTranslation, currRotation, errorArray);

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
    }

    cout<<"minMatchError "<<minMatchError<<endl;
    cout<<"final Rotation matrix "<<rotation<<endl;
    cout<<"final translation vector "<<translation<<endl;
    cout<<"bestDeltaLonLat="<<bestDeltaLonLat<<endl;


    ComputeAssembledImage(foreDEM, foreDEMGeo, backDEM, backDEMGeo,
			  assembledDEMFilename, 0, translation, rotation, bestDeltaLonLat);

  }
 
  //DRG assembler
  if (mode.compare("DRG")==0){
    Vector2 bestDeltaLonLat;
    bestDeltaLonLat(0) = 0;
    bestDeltaLonLat(1) = 0;
   
    string backDRGFilename = backFile;//"../MSLData/Mars/MER_HIRISE/PSP_001777_1650_1m_o-crop-geo.tif";
    string foreDRGFilename = foreFile;//"../MSLData/Mars/MER_HIRISE/Photo-mod.tif";
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
 
    ComputeAssembledImage(foreDRG, foreDRGGeo, backDRG, backDRGGeo,
			  assembledDRGFilename, 1, translation, rotation, bestDeltaLonLat);
  
   }
}


















