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
    string backDEMFilename = backFile;//"../MSLData/Mars/MER_HIRISE/DTEEC_001513_1655_001777_1650_U01.tif";
    string foreDEMFilename = foreFile;//"../MSLData/Mars/MER_HIRISE/Height-Sol-855.tif";

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
 
    float minMatchError = 100000000.0;
    
    /*
    //test code to determine the offset in number of background pixels for a specific lon lat offset - START
    //determine  a reasonable pixel location;
    Vector2 centerPixel;
    centerPixel(0) = backDEM.cols()/2;
    centerPixel(1) = backDEM.rows()/2;
    Vector2 origLonLat = backDEMGeo.pixel_to_lonlat(centerPixel);
    Vector2 origPixel = backDEMGeo.lonlat_to_pixel(origLonLat);
    Vector2 deltaLonLat(0.0001,0.0001);
    Vector2 newLonLat = origLonLat + deltaLonLat;
    Vector2 newPixel = backDEMGeo.lonlat_to_pixel(newLonLat);
    Vector2 deltaPixel = newPixel-origPixel;
    cout<<"deltaPixel = "<<deltaPixel<<endl;
    //test code to determine the offset in number of background pixels for a specific lon lat offset - END 
    */
    if (settings.matchingMode != 0){
       
       float matchError;
       Vector2 delta_lonlat;        
       int numVerRestarts = sqrt(settings.maxNumStarts); 
       int numHorRestarts = sqrt(settings.maxNumStarts);
      
       for (int k = -(numVerRestarts-1)/2; k < (numVerRestarts+1)/2; k++){
	 delta_lonlat(0) = k*0.0001; //~5 meters increments
	   for (int l = -(numHorRestarts-1)/2; l < (numHorRestarts+1)/2; l++){
	     delta_lonlat(1) = l*0.0001; //~5 meters increments
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
		             currTranslation, currRotation, center, matchError/*errorArray*/);
               cout<<"done."<<endl;
               /* 
	       float matchError = 0;
	       for (int m = 0; m < errorArray.size(); m++){
		 matchError = matchError+errorArray[m];
	       }
	       matchError = matchError/errorArray.size();
               cout<<"currMatchError="<<matchError<<", minMatchError="<<minMatchError<<endl;
               */

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


    vector<struct tilingParams> tileParamsArray;
    int tileSize = 128;
    int numTiles = 128;
    Vector4 paddingParams;
    paddingParams(0) = 1;//left
    paddingParams(1) = 1;//top
    paddingParams(2) = 1;//right
    paddingParams(3) = 1;//bottom
    
    Vector2 DEMOffset;
    DEMOffset(0) = (tileSize*numTiles - backDEM.cols())/2;
    DEMOffset(1) = (tileSize*numTiles - backDEM.rows())/2;
    cout<<"DEMOffset="<<DEMOffset<<endl;

    ComputeBoundariesDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo, bestDeltaLonLat, tileSize, DEMOffset, paddingParams, tileParamsArray);
    
    for (int i= 0; i <tileParamsArray.size(); i++){
        stringstream ss;
        ss<<tileParamsArray[i].horTileIndex<<"_"<<tileParamsArray[i].verTileIndex;
	string assembledDEMFilename =  resDir+"/assembled_"+ss.str()+"_dem.tif";
        cout<<assembledDEMFilename<<endl;
        
	string assembledAccFilename =  resDir+"/assembled_"+ss.str()+"_acc.tif";
        cout<<assembledDEMFilename<<endl;

	ComputeAssembledDEM(foreDEM, foreDEMGeo, backDEM, backDEMGeo, assembledDEMFilename, assembledAccFilename, 
			    translation, rotation, center, bestDeltaLonLat, tileParamsArray[i]);   
	
        //Save the assembled DEM to PointCloud file such that it can be displayed with pc_vis
	//used for debug only
	string assembledPCFilename =  resDir+"/assembled_"+ss.str()+"_pc.txt";
        cout<<assembledPCFilename<<endl;
        SaveAssembledPC(assembledDEMFilename, assembledPCFilename);
    }
   
  }

  //DRG assembler
  if (mode.compare("DRG")==0){
  
      string backDRGFilename = backFile;//"../MSLData/Mars/MER_HIRISE/DT1EA_001513_1655_001777_1650_U01.tif";
      string foreDRGFilename = foreFile;//"../MSLData/Mars/MER_HIRISE/Photo-Sol-855.tif";
      int tileSize = 512;
      int numTiles = 128;
      Vector4 paddingParams;
      paddingParams(0) = 1;//left 
      paddingParams(1) = 1;//top 
      paddingParams(2) = 2;//right
      paddingParams(3) = 2;//bottom

      Vector2 DRGOffset;
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
      
      ComputeBoundariesDRG(foreDRG, foreDRGGeo, backDRG, backDRGGeo, tileSize, DRGOffset, paddingParams, tileParamsArray);
      for (int i= 0; i <tileParamsArray.size(); i++){
        stringstream ss;
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


















