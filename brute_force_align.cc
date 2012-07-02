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

vector<vector<AlignedLOLAShot> > align_to_image(vector<vector<LOLAShot> > & trackPts, string & inputCubFile,
		vector<Matrix3x3> & trackTransforms, int transSearchWindow=0, int transSearchStep=0,
		float thetaSearchWindow=0.0, float thetaSearchStep=0.0, string image_file = "")
{
	boost::shared_ptr<DiskImageResource> rsrc(new DiskImageResourceIsis(inputCubFile));
	DiskImageView<PixelGray<float> > cub(rsrc);
	ImageView<PixelGray<float> > cubImage(cub.cols(), cub.rows());
	double nodataVal = rsrc->nodata_read();
	cubImage = apply_mask(normalize(create_mask(cub,nodataVal)),0);

	camera::IsisCameraModel model(inputCubFile);
	Vector3 center_of_moon(0,0,0);
	Vector2 pixel_location = model.point_to_pixel( center_of_moon );
	Vector3 cameraPosition = model.camera_center( pixel_location );
	Vector3 lightPosition = model.sun_position( pixel_location );
	
	//initialization step for LIMA - START	
	GetAllPtsFromCub(trackPts, model, cubImage);

	ComputeAllReflectance(trackPts, cameraPosition, lightPosition);

	transform_tracks_by_matrices(trackPts, trackTransforms);
	
	vector<vector< AlignedLOLAShot> > aligned = initialize_aligned_lola_shots(trackPts);
	Matrix3x3 trans(1, 0, 0, 0, 1, 0, 0, 0, 1);
	transform_tracks(aligned, trans, cubImage);
	float error = compute_transform_error(aligned);
	
	if (transSearchStep > 0 && thetaSearchStep > 0.0)
	{
		printf("Brute force alignment.\n");
		trans = find_tracks_transform(aligned, cubImage, 
			transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep);
	}
	for (unsigned int i = 0; i < aligned.size(); i++)
		trackTransforms[i] = gauss_newton_track(aligned[i], cubImage, trans) * trackTransforms[i];
	printf("Initial Error: %g Final Error: %g\n", error, compute_transform_error(aligned));
	
	if (image_file.length() > 0)
		SaveReflectanceImages(aligned, cubImage, image_file);

	return aligned;
}
	
int main( int argc, char *argv[] )
{
	// command line options
	string inputCSVFilename; 
	
	std::string inputCubFile;
	std::string outputFile, dataFile, imageFile, startMatricesFilename, pyramidFile;
	int transSearchWindow = 0, transSearchStep = 0;
	float thetaSearchWindow = 0.0, thetaSearchStep = 0.0;

	po::options_description general_options("Options");
	general_options.add_options()
	("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
	("inputCubFile,i", po::value<std::string>(&inputCubFile))
	("imagePyramid,p", po::value<std::string>(&pyramidFile))
	("outputFile,o", po::value<std::string>(&outputFile))
	("dataFile,d", po::value<std::string>(&dataFile))
	("outputImage", po::value<std::string>(&imageFile))
	("startMatrices,s", po::value<std::string>(&startMatricesFilename))
	("transSearchWindow", po::value<int>(&transSearchWindow))
	("transSearchStep", po::value<int>(&transSearchStep))
	("thetaSearchWindow", po::value<float>(&thetaSearchWindow))
	("thetaSearchStep", po::value<float>(&thetaSearchStep))
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

	vector<vector<LOLAShot> > trackPts = LOLAFileRead(inputCSVFilename);
	vector<Matrix3x3> trackTransforms;
	if (vm.count("startMatrices"))
		trackTransforms = load_track_transforms(startMatricesFilename);
	else
		for (unsigned int i = 0; i < trackPts.size(); i++)
			trackTransforms.push_back(Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1));
	vector<vector< AlignedLOLAShot> > aligned;

	if (vm.count("inputCubFile") > 0)
	{
		aligned = align_to_image(trackPts, inputCubFile, trackTransforms, 
				transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep, imageFile);
	}
	else if (vm.count("imagePyramid") > 0)
	{
		FILE* f = fopen(pyramidFile.c_str(), "r");
		bool first = true;
		float last_zoom;
		while (true)
		{
			float zoom_factor;
			char image_file[100];
			int ret = fscanf(f, "%g %s\n", &zoom_factor, image_file);
			if (ret != 2)
				break;
			string im(image_file);
			if (!first) // transform previous matrices
			{
				float r = last_zoom / zoom_factor;
				Matrix3x3 t1(r, 0, 0, 0, r, 0, 0, 0, 1);
				Matrix3x3 t2(1.0/r, 0, 0, 0, 1.0/r, 0, 0, 0, 1);
				for (unsigned int i = 0; i < trackTransforms.size(); i++)
					trackTransforms[i] = t1 * trackTransforms[i] * t2;
			}
			last_zoom = zoom_factor;
			if (first) // do brute force search the first time
			{
				aligned = align_to_image(trackPts, im, trackTransforms, 
					transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep);
				first = false;
			}
			else if (zoom_factor == 1) // last level, save image
			{
				aligned = align_to_image(trackPts, im, trackTransforms, 
					0.0, 0.0, 0.0, 0.0, imageFile);
			}
			else
			{
				aligned = align_to_image(trackPts, im, trackTransforms);
			}
		}
		fclose(f);
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
	for (unsigned int i = 0; i < trackTransforms.size(); i++)
	{
		Matrix3x3 & m = trackTransforms[i];
		fprintf(output, "%g %g %g %g %g %g\n", m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2));
	}
	fclose(output);

	if (dataFile.length() > 0)
		save_track_data(aligned, dataFile);
 
	return 0;

}

