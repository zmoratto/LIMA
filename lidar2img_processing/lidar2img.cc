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
#include "../util.h"
#include "../tracks.h"
#include "../match.h"
#include "../coregister.h"
#include "../display.h"
#include "../weights.h"
#include "../featuresLOLA.h"

// just initialize the tracks, don't align them
vector<vector<AlignedLOLAShot> > setup_tracks(vector<vector<LOLAShot> > & trackPts, const string & inputCubFile,
		Matrix3x3 trans)
{
	DiskImageResourceIsis rsrc(inputCubFile);
	ImageView<PixelGray<float> > cubImage;
	read_image(cubImage, rsrc);
	//double nodataVal = rsrc.nodata_read();
	cubImage = normalize(cubImage);//apply_mask(normalize(create_mask(cubImage,nodataVal)),0);

	camera::IsisCameraModel model(inputCubFile);
	Vector3 center_of_moon(0,0,0);
	Vector2 pixel_location = model.point_to_pixel( center_of_moon );
	Vector3 cameraPosition = model.camera_center( pixel_location );
	Vector3 lightPosition = model.sun_position( pixel_location );
	
	//initialization step for LIMA - START	
	GetAllPtsFromCub(trackPts, model, cubImage);

	ComputeAllReflectance(trackPts, cameraPosition, lightPosition);

	vector<vector< AlignedLOLAShot> > aligned = initialize_aligned_lola_shots(trackPts);
	transform_tracks(aligned, trans, cubImage);
	return aligned;
}

vector<vector<AlignedLOLAShot> > align_to_image(vector<vector<LOLAShot> > & trackPts, ImageView<PixelGray<float> > cubImage,
		vector<Matrix3x3> & trackTransforms, bool globalAlignment=true, int transSearchWindow=0, int transSearchStep=0,
		float thetaSearchWindow=0.0, float thetaSearchStep=0.0, string image_file = "")
{
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
	else if (globalAlignment)
	{
		printf("Global Gauss-Newton alignment.\n");
		vector<AlignedLOLAShot> oneTrack;
		for (unsigned int i = 0; i < aligned.size(); i++)
			oneTrack.insert(oneTrack.end(), aligned[i].begin(), aligned[i].end());
		trans = gauss_newton_track(oneTrack, cubImage, trackTransforms[0]);
		transform_tracks(aligned, trans, cubImage);
		for (unsigned int i = 0; i < trackTransforms.size(); i++)
			trackTransforms[i] = trans;// * trackTransforms[i];
	}
	else
	{
		for (unsigned int i = 0; i < aligned.size(); i++)
			trackTransforms[i] = gauss_newton_track(aligned[i], cubImage, trans) * trackTransforms[i];
	}
	printf("Initial Error: %g Final Error: %g\n", error, compute_transform_error(aligned));
	
	if (image_file.length() > 0)
		SaveReflectanceImages(aligned, cubImage, image_file);

	return aligned;
}

vector<vector<AlignedLOLAShot> > align_to_image_pyramid(vector<vector<LOLAShot> > & trackPts, const string & image_file,
		vector<Matrix3x3> & trackTransforms, string outputImage = "")
{
	int ZOOM_MULTIPLIER = 2;
	
	int zoom_factor = 8;

	vector<vector< AlignedLOLAShot> > aligned;

	bool first = true;
	DiskImageResourceIsis rsrc(image_file);
	ImageView<PixelGray<float> > cubImage;
	read_image(cubImage, rsrc);
	//double nodataVal = rsrc.nodata_read();
	cubImage = normalize(cubImage);//apply_mask(normalize(create_mask(cubImage,nodataVal)),0);
	
	camera::IsisCameraModel model(image_file);
	Vector3 center_of_moon(0,0,0);
	Vector2 pixel_location = model.point_to_pixel( center_of_moon );
	Vector3 cameraPosition = model.camera_center( pixel_location );
	Vector3 lightPosition = model.sun_position( pixel_location );
	
	GetAllPtsFromCub(trackPts, model, cubImage);
	ComputeAllReflectance(trackPts, cameraPosition, lightPosition);

	transform_tracks_by_matrix(trackPts, Matrix3x3(1.0 / zoom_factor, 0, 0, 0, 1.0 / zoom_factor, 0, 0, 0, 0));

	while (zoom_factor >= 1.0)
	{
		if (!first) // transform previous matrices
		{
			Matrix3x3 t1((float)ZOOM_MULTIPLIER, 0, 0, 0, (float)ZOOM_MULTIPLIER, 0, 0, 0, 1);
			Matrix3x3 t2(1.0/ZOOM_MULTIPLIER, 0, 0, 0, 1.0/ZOOM_MULTIPLIER, 0, 0, 0, 1);
			for (unsigned int i = 0; i < trackTransforms.size(); i++)
				trackTransforms[i] = t1 * trackTransforms[i] * t2;
			transform_tracks_by_matrix(trackPts, Matrix3x3(ZOOM_MULTIPLIER, 0, 0, 0, ZOOM_MULTIPLIER, 0, 0, 0, 1));
		}
		else
		{
			first = false;
		}
		if (zoom_factor == 1) // last level, save image
		{
			aligned = align_to_image(trackPts, cubImage, trackTransforms, true,
				0.0, 0.0, 0.0, 0.0, outputImage);
		}
		else
		{
			ImageView<PixelGray<float> > img = resize(cubImage, cubImage.rows() / zoom_factor, cubImage.cols() / zoom_factor);
			
			aligned = align_to_image(trackPts, img, trackTransforms);
		}
		zoom_factor /= ZOOM_MULTIPLIER;
	}

	return aligned;
}

int main( int argc, char *argv[] )
{
	// command line options
	string inputCSVFilename; 
	
	std::string inputCubFile, tracksListFile, gcpFile;
	std::string outputFile, dataFile, imageFile;

	po::options_description general_options("Options");
	general_options.add_options()
	("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
	("tracksList,t", po::value<std::string>(&tracksListFile))
	("inputCubFile,i", po::value<std::string>(&inputCubFile))
	("outputFile,o", po::value<std::string>(&outputFile))
	("dataFile,d", po::value<std::string>(&dataFile))
	("outputImage", po::value<std::string>(&imageFile))
	("gcpFile,g", po::value<std::string>(&gcpFile))
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
	if (vm.count("Lidar-filename"))
		trackPts = LOLAFileRead(inputCSVFilename);
	else if (vm.count("tracksList"))
		trackPts = LOLAFilesRead(tracksListFile);
	else
	{
		fprintf(stderr, "No track files specified.\n");
		return 1;
	}
	vector<Matrix3x3> trackTransforms;
	for (unsigned int i = 0; i < trackPts.size(); i++)
		trackTransforms.push_back(Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1));
	
	// use this later, only if saving to file
	vector<gcp> gcpArray;
	if (gcpFile.length() > 0)
		gcpArray = ComputeSalientLOLAFeatures(trackPts);

	vector<vector< AlignedLOLAShot> > aligned;

	if (vm.count("inputCubFile") > 0)
	{
		aligned = align_to_image_pyramid(trackPts, inputCubFile, trackTransforms, imageFile);
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
		for (unsigned int i = 0; i < trackTransforms.size(); i++)
		{
			Matrix3x3 & m = trackTransforms[i];
			fprintf(output, "%g %g %g %g %g %g\n", m(0, 0), m(0, 1), m(0, 2), m(1, 0), m(1, 1), m(1, 2));
		}
		fclose(output);
	}
	
	if (gcpFile.length() > 0)
	{
		UpdateGCP(aligned, inputCubFile, gcpArray);
		SaveGCPoints(gcpArray, gcpFile);
	}

	if (dataFile.length() > 0)
		save_track_data(aligned, dataFile);
 
	return 0;

}

