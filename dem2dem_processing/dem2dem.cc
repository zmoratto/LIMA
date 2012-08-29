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
#include "icp.h"

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
  if( ReadConfigFile(configFilename, &settings) )
	{
	std::cerr << "Config file " << configFilename << " found." << endl;
	}
  else
	{
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
  
  string backDEMFilename = backFile;
  string foreDEMFilename = foreFile;
   
  //open the ref or model DEM 
  cout<<"opening"<<backDEMFilename<<"..."<<endl; 
  //DiskImageView<PixelGray<float> >  backDEM(backDEMFilename);

  boost::shared_ptr<DiskImageResource> back_rsrc( new DiskImageResourceGDAL(backDEMFilename) );
  if (back_rsrc->has_nodata_read()){
      settings.noDataVal = back_rsrc->nodata_read();
      //if( verbose > 0 ){
	 cout << "Found nodata value, " << settings.noDataVal << ", in " << backDEMFilename << endl;
	 //}
  }
  else /*if( verbose > 0 )*/{
	 cout << "Using default nodata value: " << settings.noDataVal << endl;
  }
  DiskImageView<PixelGray<float> > backDEM( back_rsrc );



  GeoReference backDEMGeo;
  read_georeference(backDEMGeo, backDEMFilename);
  float back_radius = backDEMGeo.datum().semi_major_axis();
  cout<<"back_radius="<<back_radius<<endl;
  cout<<"done"<<endl;
  
  //open the matching DEM
 
  cout<<"opening"<<foreDEMFilename<<"..."<<endl;  
  //DiskImageView<PixelGray<float> >  foreDEM(foreDEMFilename);
 
  boost::shared_ptr<DiskImageResource> fore_rsrc( new DiskImageResourceGDAL(foreDEMFilename) );
  if (fore_rsrc->has_nodata_read()){
      settings.noDataVal = fore_rsrc->nodata_read();
      //if( verbose > 0 ){
	 cout << "Found nodata value, " << settings.noDataVal << ", in " << foreDEMFilename << endl;
	 //}
  }
  else /*if( verbose > 0 )*/{
	 cout << "Using default nodata value: " << settings.noDataVal << endl;
  }
  DiskImageView<PixelGray<float> > foreDEM( fore_rsrc );
 
  
  GeoReference foreDEMGeo;
  read_georeference(foreDEMGeo, foreDEMFilename);
  float fore_radius = foreDEMGeo.datum().semi_major_axis();
  cout<<"fore_radius="<<fore_radius<<endl;
  cout<<"done"<<endl;
   
  float minMatchError = 100000000.0;
     
  Vector2 delta_lonlat; 
       
  printf("feature extraction ...\n");
	    
  vector<Vector3> featureArray = GetFeatures(foreDEM, foreDEMGeo, backDEM, backDEMGeo, 
					     settings.samplingStep, delta_lonlat, settings.noDataVal);
  vector<float> errorArray;
  errorArray.resize(featureArray.size());
  
  Vector3 currTranslation;
  Matrix<float, 3,3 > currRotation;
  ICP_DEM_2_DEM(featureArray, backDEM, backDEMGeo, foreDEMGeo, settings,
		currTranslation, currRotation, center, errorArray);

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
	       
  cout<<"minMatchError "<<minMatchError<<endl;
  cout<<"final Rotation matrix "<<rotation<<endl;
  cout<<"final translation vector "<<translation<<endl;
  cout<<"bestDeltaLonLat="<<bestDeltaLonLat<<endl;  


  //TO DO: compute foreDEM
  ImageView<PixelGray<float> > transfForeDEM(foreDEM.cols(), foreDEM.rows());
  for (unsigned int i=0; i < foreDEM.cols(); i++){
    for (unsigned int j=0; j < foreDEM.rows(); j++){
	transfForeDEM(i,j)=-10000; 
      }
  }
  
  for (unsigned int i=0; i < foreDEM.cols(); i++){
    for (unsigned int j=0; j < foreDEM.rows(); j++){
      //start
         Vector2 forePix(i,j);
         Vector2 fore_lon_lat = foreDEMGeo.pixel_to_lonlat(forePix);
      
         //spherical to cartesian
	 Vector3 fore_lon_lat_rad;
	 fore_lon_lat_rad(0) = fore_lon_lat(0);
	 fore_lon_lat_rad(1) = fore_lon_lat(1);
	 fore_lon_lat_rad(2) = foreDEM(i,j);
	 Vector3 fore_xyz = foreDEMGeo.datum().geodetic_to_cartesian(fore_lon_lat_rad); 
	 
	 //Vector3 transf_fore_xyz = rotation*(fore_xyz-foreCenter_xyz)+translation+foreCenter_xyz;
         Vector3 transf_fore_xyz = fore_xyz + translation;

	 //transform back in spherical coords
	 Vector3 transf_fore_lon_lat_rad = foreDEMGeo.datum().cartesian_to_geodetic(transf_fore_xyz);
	 Vector2 transf_fore_lon_lat;
	 transf_fore_lon_lat(0)=transf_fore_lon_lat_rad(0);
	 transf_fore_lon_lat(1)=transf_fore_lon_lat_rad(1);
        
	 Vector2 backPix= foreDEMGeo.lonlat_to_pixel(transf_fore_lon_lat);
	 int x = backPix(0); 
         int y = backPix(1);
         
	 //overwrite the assembled image 
	 if ((x>0) && (y>0) && (x<transfForeDEM.cols()) && (y<transfForeDEM.rows())){ 
	    transfForeDEM(x,y) =  transf_fore_lon_lat_rad(2);
	 }
   
    //end 
      }
  }
  
  //save fore DEM to file.  
  string transfForeDEMFilename =  resDir+"/transf_dem.tif";
  cout<<transfForeDEMFilename<<endl;
  write_georeferenced_image(transfForeDEMFilename,
                            transfForeDEM,
                            foreDEMGeo, TerminalProgressCallback("{Core}","Processing:"));
  
}


















