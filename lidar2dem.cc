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
  std::string configFilename="lidar2dem_settings.txt";
  std::vector<std::string> DEMFiles;
  std::vector<std::string> cubFiles;
  std::string resDir = "../results";
 
  po::options_description general_options("Options");
  general_options.add_options()
    ("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
    ("DEMFiles,d", po::value<std::vector<std::string> >(&DEMFiles))
    ("results-directory,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
    ("settings-filename,s", po::value<std::string>(&configFilename)->default_value("lidar2dem_settings.txt"), "settings filename.")
    ("help,h", "Display this help message");
  

  po::options_description hidden_options("");

  hidden_options.add_options()
    ("cubFiles,c", po::value<std::vector<std::string> >(&DEMFiles));
 
  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("DEMFiles", -1);

  std::ostringstream usage;
  usage << "Description: main code for Lidar to DEM co-registration" << std::endl << std::endl;
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

  if(( vm.count("DEMFiles") < 1 )) {
    std::cerr << "Error: Must specify at least  one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }
 
  struct CoregistrationParams settings;
  ReadConfigFile((char*)configFilename.c_str(), &settings);
  PrintGlobalParams(&settings);

  //read LOLA tracks
  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);


  for (int index = 0; index <  DEMFiles.size(); index++){

      string inputDEMFilename = DEMFiles[index];
      cout<<"inputDEMFilename"<<inputDEMFilename<<endl;

      //create the results directory and prepare the output filenames - START
      system("mkdir ../results");

      string DEMFilenameNoPath = sufix_from_filename(inputDEMFilename);
      
      string altitudePtsFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "alt.txt";
      string demPtsFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "dem.txt";
      string lolaTracksFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "_lola.tif";  
      string lolaFeaturesFilename = "../results" + prefix_less3_from_filename(DEMFilenameNoPath) + "_features_lola.txt";  
      //create the results directory and prepare the output filenames - END

  
      //read DEM file
      
      boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(inputDEMFilename) );
      double nodata_value;
      if (rsrc->has_nodata_read()){
	nodata_value = rsrc->nodata_read();
      }
      else{
	nodata_value = settings.noDataVal;
      }
      cout<<"nodata value="<<nodata_value<<endl;
 
      DiskImageView<PixelGray<float> > DEM( rsrc );
     

      GeoReference DEMGeo;
      read_georeference(DEMGeo, inputDEMFilename);
 
      ImageViewRef<PixelGray<float> >   interpDEM = interpolate(edge_extend(DEM.impl(),
									    ConstantEdgeExtension()),
								BilinearInterpolation());

      //select DEM points closes to LOLA tracks
      GetAllPtsFromDEM(trackPts, interpDEM, DEMGeo, nodata_value);

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
      vector<Vector3> featureArray;//DEM
      vector<Vector3> modelArray;//LOLA
      vector<float> errorArray;//error array same size as model and feature
   

      //copy info to featureArray and modelArray
      for(int k = 0; k < trackPts.size();k++){
	for(int i = 0; i < trackPts[k].size(); i=i+settings.samplingStep(0)){
       
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

      //run ICP-matching
      ICP(featureArray, modelArray, settings,currTranslation, currRotation, errorArray);
      
      cout<<"Translation="<<currTranslation<<endl;
      cout<<"Rotation="<<currRotation<<endl;
  }

  return 0;
}


















