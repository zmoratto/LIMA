/**
 * This program takes as input a cub file which is an image of a planetary surface,
 * and a CSV file of LOLA lidar points. It performs a brute force search over a
 * local area to find the transform which best fits the lidar points to the image
 * according to their estimate reflectance.
 **/

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
#include <vw/Mosaic/ImageComposite.h>
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
#include "lidar_tracks/tracks.h"
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"
#include "featuresLOLA.h"

int main( int argc, char *argv[] )
{
	// command line options
	string inputCSVFilename; 
	
	std::string inputCubFile, tracksListFile, gcpFile;
	std::string outputFile, dataFile, imageFile;

	po::options_description general_options("Options");
	general_options.add_options()
	("lidarFile,l", po::value<std::string>(&inputCSVFilename))
	("tracksList,t", po::value<std::string>(&tracksListFile))
	("inputCubFile,i", po::value<std::string>(&inputCubFile))
	("outputFile,o", po::value<std::string>(&outputFile))
	("dataFile,d", po::value<std::string>(&dataFile))
	("outputImage", po::value<std::string>(&imageFile))
	("gcpDir,g", po::value<std::string>(&gcpFile))
	("help,h", "Display this help message");
	
	po::options_description options("Allowed Options");
	options.add(general_options);

	po::positional_options_description p;
	p.add("inputCubFile", -1);

	std::ostringstream usage;
	usage << "Description: main code for Lidar to image co-registration" << std::endl << std::endl;
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

	vector<vector<LOLAShot> > trackPts;
	if (vm.count("lidarFile"))
		trackPts = LOLAFileRead(inputCSVFilename);
	else if (vm.count("tracksList"))
		trackPts = LOLAFilesRead(tracksListFile);
	else
	{
		fprintf(stderr, "No track files specified.\n");
		return 1;
	}
	Matrix3x3 trans(1, 0, 0, 0, 1, 0, 0, 0, 1);
	
	// use this later, only if saving to file
	vector<gcp> gcpArray;
	if (gcpFile.length() > 0)
		gcpArray = ComputeSalientLOLAFeatures(trackPts);

	vector<vector< AlignedLOLAShot> > aligned;

	if (vm.count("inputCubFile") > 0)
	{
		aligned = align_to_image_pyramid(trackPts, inputCubFile, trans, imageFile);
	}
	else
	{
		fprintf(stderr, "Must specify either a cub file or an image pyramid file!\n");
		return 1;
	}
	
	FILE* output = stdout;
	if (outputFile.length() > 0)
	{
		output = fopen(outputFile.c_str(), "w");
		if (output == NULL)
		{
			fprintf(stderr, "Failed to open output file %s.\n", outputFile.c_str());
			output = stdout;
		}
	}
	fprintf(output, "%g %g %g %g %g %g\n", trans(0, 0), trans(0, 1), trans(0, 2), trans(1, 0), trans(1, 1), trans(1, 2));
	if (outputFile.length() > 0)
		fclose(output);
	
	if (gcpFile.length() > 0)
	{
		UpdateGCP(aligned, inputCubFile, gcpArray);
		SaveGCPoints(gcpArray, gcpFile);
	}

	if (dataFile.length() > 0)
		save_track_data(aligned, dataFile);
 
	return 0;

}

