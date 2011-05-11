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
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"
#include "icp.h"

int main( int argc, char *argv[] )
{
string inputCSVFilename; 
int verbose;
std::string configFilename;
std::vector<std::string> DEMFiles;
std::vector<std::string> cubFiles;
std::string resDir;
 
po::options_description general_options("Options");
general_options.add_options()
	("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
    ("DEMFiles,d", po::value<std::vector<std::string> >(&DEMFiles))
    ("results-directory,r", 
		po::value<std::string>(&resDir)->default_value("../results"), 
		"results directory.") // Currently a no-op until we implement it below.
    ("settings-filename,s", 
		po::value<std::string>(&configFilename)->default_value("lidar2dem_settings.txt"), 
		"settings filename.")
    ("verbose,v", 
		po::value<int>(&verbose)->default_value(1), 
		"Verbosity level, zero emits no messages.")
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
try
	{
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
	}
catch ( po::error &e )
	{
    std::cout << "An error occured while parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
	}

if( vm.count("help") )
	{
    std::cerr << usage.str() << std::endl;
    return 1;
	}

if(( vm.count("DEMFiles") < 1 )) {
    std::cerr << "Error: Must specify at least  one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
}

// Set up VW logging
vw_log().console_log().rule_set().add_rule(vw::InfoMessage,"*");

 
struct CoregistrationParams settings;
if( ReadConfigFile(configFilename, &settings) && verbose > 0 )
	{
	cerr << "Config file " << configFilename << " found." << endl;
	}
else if( verbose > 0 )
	{
	cerr << "Config file " << configFilename << " not found, using defaults." << endl;
	}
//PrintGlobalParams(&settings);
if( verbose > 0 ){ cerr << settings << endl; }

  //read LOLA tracks
  vector<vector<LOLAShot> > trackPts;
  try
	{
	trackPts =  CSVFileRead(inputCSVFilename);
	if( verbose > 0 ){ cerr << "Tracks file: " << inputCSVFilename << endl; }
	}
  catch (const vw::IOErr& error)
	{
	cerr << error.what() << endl;
	exit(1);
	}


  for (unsigned int index = 0; index <  DEMFiles.size(); index++){

      string inputDEMFilename = DEMFiles[index];
      if( verbose > 0 ){ cout <<"DEM filename: " << inputDEMFilename << endl; }

	/* Re-implement this once we figure out what needs to be written out.
      //create the results directory and prepare the output filenames - START
      system( resDir.insert(0,"mkdir ").c_str() );

      string DEMFilenameNoPath = sufix_from_filename(inputDEMFilename);
      
      string altitudePtsFilename = resDir + prefix_less3_from_filename(DEMFilenameNoPath) + "alt.txt";
      string demPtsFilename = resDir + prefix_less3_from_filename(DEMFilenameNoPath) + "dem.txt";
      string lolaTracksFilename = resDir + prefix_less3_from_filename(DEMFilenameNoPath) + "_lola.tif";  
      string lolaFeaturesFilename = resDir + prefix_less3_from_filename(DEMFilenameNoPath) + "_features_lola.txt";  
      //create the results directory and prepare the output filenames - END
	*/ 

  
      //read DEM file
      
      boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(inputDEMFilename) );
      if (rsrc->has_nodata_read()){
        settings.noDataVal = rsrc->nodata_read();
		if( verbose > 0 )
			{
			cout 	<< "Found nodata value, " << settings.noDataVal 
					<< ", in " << inputDEMFilename << endl;
			}
      }
      else if( verbose > 0 ){
		cout << "Using default nodata value: " << settings.noDataVal << endl;
      }
 
      DiskImageView<PixelGray<float> > DEM( rsrc );
     
      GeoReference DEMGeo;
      read_georeference(DEMGeo, inputDEMFilename);
 
      ImageViewRef<PixelGray<float> >   interpDEM = interpolate(edge_extend(DEM.impl(),
									    ConstantEdgeExtension()),
								BilinearInterpolation());

      //select DEM points closes to LOLA tracks
      GetAllPtsFromDEM(trackPts, interpDEM, DEMGeo, settings.noDataVal);

      //initialization step for LIDEM - END
      /*
      if (settings.analyseFlag == 1){
         SaveDEMPoints(trackPts, inputDEMFilename, demPtsFilename);
         int numVerPts = 6000;
         int numHorPts = 6000;
         MakeGrid(trackPts, numVerPts, numHorPts, lolaTracksFilename, trackIndices);
      }
      */
  
      Vector3 currTranslation;
      Matrix<float, 3,3 > currRotation;
      vector<Vector3> featureArray;//DEM
      vector<Vector3> modelArray;//LOLA
      vector<float> errorArray;//error array same size as model and feature
      currRotation[0][0] = 1.0;
      currRotation[1][1] = 1.0;
      currRotation[2][2] = 1.0;
      Vector3 center;

      //copy info to featureArray and modelArray
      for(unsigned int k = 0; k < trackPts.size();k++){
	for(unsigned int i = 0; i < trackPts[k].size(); i=i+settings.samplingStep(0)){
       
	  if ((trackPts[k][i].valid == 1) && (trackPts[k][i].DEMPt[2].valid==1)){
	    Vector3 model;
	    // Vector3 feature;
         
	    //this is the LOLA data
        //    float radius = DEMGeo.datum().semi_major_axis();
	    model = trackPts[k][i].LOLAPt[2].coords;

		model.z() *= 1000; // LOLA data is in km, DEMGeo is in m (for LROC DTMs).
          

	    if ((model[0] >1e-100) && (model[1] >1e-100) && (model[2]>1e-100) ){
			if( verbose > 1 ){ cout<<"altitude="<<model[2]<<endl; }
	      /*
              feature[0] = trackPts[k][i].LOLAPt[2].coords(0); 
	      feature[1] = trackPts[k][i].LOLAPt[2].coords(1); 
	      feature[2] = trackPts[k][i].DEMPt[2].val;
	      */
	      //copy the model and features into cartesian coordinates
	      //feature = DEMGeo.datum().geodetic_to_cartesian(feature);

	      model = DEMGeo.datum().geodetic_to_cartesian(model);
         	 
	      //featureArray.push_back(feature);
	      modelArray.push_back(model);
	    }
	  }
	}
      }

      featureArray.resize(modelArray.size());
      if( verbose > 0 )
		{
		cout << "Array Sizes: " << modelArray.size() <<" and "<< featureArray.size() << endl;
		}

      //run ICP-matching
      ICP_LIDAR_2_DEM(featureArray, interpDEM, DEMGeo, modelArray, settings, 
                      currTranslation, currRotation, center, errorArray);
     
      if( verbose >= 0 )
		{ 
			cout<<"Translation="<<currTranslation<<endl;
			cout<<"Rotation="<<currRotation<<endl;
			cout<<"Center="<<center<<endl;
		}
  }

  return 0;
}
