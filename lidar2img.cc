// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::camera;
using namespace std;

#include <math.h>
#include "util.h"
#include "tracks.h"
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"
#include "featuresLOLA.h"
	
std::vector<int> determine_overlapping(vector<vector<LOLAShot> > trackPts, std::vector<std::string> camCubFiles)
{
	std::vector<int> overlapIndices;
	Vector4 lat_lon_bb = FindMinMaxLat(trackPts); 
	Vector4 lon_lat_bb;
	lon_lat_bb[0]=lat_lon_bb[2];
	lon_lat_bb[1]=lat_lon_bb[3];
	lon_lat_bb[2]=lat_lon_bb[0];
	lon_lat_bb[3]=lat_lon_bb[1];
	overlapIndices = makeOverlapList(camCubFiles, lon_lat_bb);

	return overlapIndices;
}

int main( int argc, char *argv[] )
{
	// command line options
	string inputCSVFilename; 
	std::string configFilename="lidar2img_settings.txt";
	
	std::vector<std::string> camCubFiles;
	std::string resDir = "../results";
	std::string mapCubDir = "../data/map";
	std::string drgDir = "../data/drg";
	std::vector<std::string> inputCubFiles;
	struct CoregistrationParams settings;

	po::options_description general_options("Options");
	general_options.add_options()
	("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
	("inputCubFiles,i", po::value<std::vector<std::string> >(&inputCubFiles))
	("mapCub-directory,m", po::value<std::string>(&mapCubDir)->default_value("../data/map"), "map cub directory.")
	("drg-directory,d", po::value<std::string>(&drgDir)->default_value("../data/drg"), "drg directory.")
	("results-directory,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
	("settings-filename,s", po::value<std::string>(&configFilename)->default_value("lidar2img_settings.txt"), "settings filename.")
	("help,h", "Display this help message");
	
	po::options_description options("Allowed Options");
	options.add(general_options);

	po::positional_options_description p;
	//p.add("camCubFileList", -1);
	p.add("inputCubFiles", -1);

	std::ostringstream usage;
	usage << "Description: main code for Lidar to image co-registration" << std::endl << std::endl;
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

	if(( vm.count("inputCubFiles") < 1 )) {
	std::cerr << "Error: Must specify at least one cub image file!" << std::endl << std::endl;
	std::cerr << usage.str();
	return 1;
	}
	
	//#if 0
	if( ReadConfigFile(configFilename, &settings) )
		std::cerr << "Config file " << configFilename << " found." << endl;
	else
		std::cerr << "Config file " << configFilename << " not found, using defaults." << endl;
	//PrintGlobalParams(&settings);
	//std::cerr << settings << endl;

	
	//determine if inputCubFiles is a text file containing a list of .cub files, one .cub file or a set of .cub files
	//by reading the file extension. A text file containing the DEM list *must* have extension .txt 
	camCubFiles = AccessDataFilesFromInput(inputCubFiles);
	

	//create the results directory and prepare the output filenames - START

	string makeResDirCmd = "mkdir -p " + resDir;
	int ret = system(makeResDirCmd.c_str()); 
	if (ret) exit(1);

	vector<vector<LOLAShot> > trackPts =	CSVFileRead(inputCSVFilename);
	std::vector<int> overlapIndices = determine_overlapping(trackPts, camCubFiles);
	vector<gcp> gcpArray = ComputeSalientLOLAFeatures(trackPts);
	
	//Save3DImage(trackPts, resDir + "/3d_moon.obj");
	
	//save the GCP
	string gcpFilenameRoot = resDir+"/"+GetFilenameNoExt(GetFilenameNoPath(inputCSVFilename));
 
	for (unsigned int k = 0; k < overlapIndices.size(); k++)
	{
		string overlapCamCubFile = camCubFiles[overlapIndices[k]];
		string overlapMapCubFile = GetFilenameNoPath(overlapCamCubFile); 
		overlapMapCubFile = mapCubDir+string("/")+overlapMapCubFile;
		FindAndReplace(overlapMapCubFile, ".cub", "_map.cub"); 
		string overlapDRGFilename = drgDir+string("/")+GetFilenameNoPath(overlapCamCubFile);
		FindAndReplace(overlapDRGFilename, ".cub", "_drg.tif");
	
		boost::shared_ptr<DiskImageResource> rsrc(new DiskImageResourceIsis(overlapCamCubFile));
		DiskImageView<PixelGray<float> > cub(rsrc);
		ImageView<PixelGray<float> > cubImage(cub.cols(), cub.rows());
		double nodataVal = rsrc->nodata_read();
		cubImage = apply_mask(normalize(create_mask(cub,nodataVal)),0);
	
	
		camera::IsisCameraModel model(overlapCamCubFile);
		//camera::IsisCameraModel model(overlapMapCubFile); // for DRG
		Vector3 center_of_moon(0,0,0);
		Vector2 pixel_location = model.point_to_pixel( center_of_moon );
		Vector3 cameraPosition = model.camera_center( pixel_location );
		Vector3 lightPosition = model.sun_position( pixel_location );
	
		//initialization step for LIMA - START	
		GetAllPtsFromCub(trackPts, model, cubImage);
		/*boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(overlapDRGFilename) );
		DiskImageView<PixelGray<uint8> > DRG( rsrc );
		GeoReference DRGGeo;
		read_georeference(DRGGeo, overlapDRGFilename);
		GetAllPtsFromImage(trackPts, DRG, DRGGeo);*/
		
		int	numValidReflPts = ComputeAllReflectance(trackPts, cameraPosition, lightPosition);
		vector<vector< AlignedLOLAShot> > aligned = initialize_aligned_lola_shots(trackPts);

		transform_tracks(aligned, Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1), cubImage);
		//initialization step for LIMA - END 
	
		
		if (numValidReflPts < 100)
			fprintf(stderr, "Not enough reflectance points, aborting.\n");

		//find_track_transforms(aligned, overlapCamCubFile);
		//Matrix3x3 trans = find_tracks_transform(aligned, overlapCamCubFile);
		Matrix3x3 trans(1, 0, 10, 0, 1, -15, 0, 0, 1);
		transform_tracks(aligned, trans, cubImage);
		printf("Best transform:\n");
		for (int i = 0; i < 3; i++)
			printf("%g %g %g\n", trans(i, 0), trans(i, 1), trans(i, 2));
		//gauss_newton_track(aligned[2], overlapCamCubFile, trans);
		//return 0;
		for (unsigned int i = 0; i < aligned.size(); i++)
			gauss_newton_track(aligned[i], cubImage, trans);
		//vector<float> initMatchingErrorArray;
		//vector<Vector4> matchArray = FindMatches2D(trackPts, overlapCamCubFile, settings.matchWindowHalfSize, 80, initMatchingErrorArray);
		std::stringstream out2;
		out2 << resDir << "/reflectance_" << k << ".tif";
		SaveReflectanceImages(aligned, cubImage, out2.str());
		std::stringstream out4;
		out4 << resDir << "/track_data_" << k << ".txt";
		save_track_data(aligned, out4.str());
		//SaveReflectanceImages(aligned, overlapDRGFilename, out2.str(), false);
	
		//save to GCP structure. 
		//UpdateGCP(trackPts, matchArray, overlapCamCubFile, gcpArray);
		
		std::stringstream out3;
		out3 << gcpFilenameRoot << "_match_" << k << ".tif";
		//SaveAdjustedReflectanceImages(gcpArray, aligned, overlapCamCubFile, out3.str(), k, false);
	}	
 
	//SaveGCPoints(gcpArray,	gcpFilenameRoot);
 
	return 0;

}

