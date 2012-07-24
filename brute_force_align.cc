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
#include "tracks.h"
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"
#include "featuresLOLA.h"

vector<vector<AlignedLOLAShot> > align_to_image(vector<vector<LOLAShot> > & trackPts, string & inputCubFile,
		vector<Matrix3x3> & trackTransforms, bool globalAlignment=false, int transSearchWindow=0, int transSearchStep=0,
		float thetaSearchWindow=0.0, float thetaSearchStep=0.0, string image_file = "")
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
	else if (globalAlignment)
	{
		printf("Global Gauss-Newton alignment.\n");
		vector<AlignedLOLAShot> oneTrack;
		for (unsigned int i = 0; i < aligned.size(); i++)
			oneTrack.insert(oneTrack.end(), aligned[i].begin(), aligned[i].end());
		trans = gauss_newton_track(oneTrack, cubImage, trans);
		transform_tracks(aligned, trans, cubImage);
		for (unsigned int i = 0; i < trackTransforms.size(); i++)
			trackTransforms[i] = trans * trackTransforms[i];
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

vector<vector<AlignedLOLAShot> > align_to_image_pyramid(vector<vector<LOLAShot> > & trackPts, const string & pyramidFile,
		vector<Matrix3x3> & trackTransforms, bool globalAlignment=false, int transSearchWindow=0, int transSearchStep=0,
		float thetaSearchWindow=0.0, float thetaSearchStep=0.0, string outputImage = "")
{
	vector<vector< AlignedLOLAShot> > aligned;

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
			if (!globalAlignment)
				aligned = align_to_image(trackPts, im, trackTransforms, true,
					transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep);
			aligned = align_to_image(trackPts, im, trackTransforms, globalAlignment,
				transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep);
			first = false;
		}
		else if (zoom_factor == 1) // last level, save image
		{
			aligned = align_to_image(trackPts, im, trackTransforms, globalAlignment,
				0.0, 0.0, 0.0, 0.0, outputImage);
		}
		else
		{
			aligned = align_to_image(trackPts, im, trackTransforms, globalAlignment);
		}
	}
	fclose(f);

	return aligned;
}

// returns a transformation matrix from coordinates in image 2 to coordinates in image 1
Matrix3x3 find_correspondence(vector<vector<AlignedLOLAShot> > aligned, vector<vector<AlignedLOLAShot> > aligned2)
{
	int num_points = 0;
	for (unsigned int i = 0; i < aligned.size(); i++)
	{
		for (unsigned int j = 0; j < aligned[i].size(); j++)
			if (aligned[i][j].image_x >= 0 && aligned2[i][j].image_x >= 0)
				num_points++;
	}
	// [a b c; d e f; 0 0 1] * [x;y;1] = [x';y';1] = [ax + by + c; dx + ey + f; 1]
	// [x y 1 0 0 0; 0 0 0 x y 1] * [a;b;c;d;e;f]
	Matrix<double> A(2 * num_points, 6);
	Matrix<double> b(2 * num_points, 1);
	int index = 0;
	for (unsigned int i = 0; i < aligned.size(); i++)
		for (unsigned int j = 0; j < aligned[i].size(); j++)
			if (aligned[i][j].image_x >= 0 && aligned2[i][j].image_x >= 0)
			{
				A(2 * index    , 0) = aligned2[i][j].image_x;
				A(2 * index    , 1) = aligned2[i][j].image_y;
				A(2 * index    , 2) = 1;
				A(2 * index    , 3) = 0;
				A(2 * index    , 4) = 0;
				A(2 * index    , 5) = 0;
				A(2 * index + 1, 0) = 0;
				A(2 * index + 1, 1) = 0;
				A(2 * index + 1, 2) = 0;
				A(2 * index + 1, 3) = aligned2[i][j].image_x;
				A(2 * index + 1, 4) = aligned2[i][j].image_y;
				A(2 * index + 1, 5) = 1;
				b(2 * index    , 0) = aligned[i][j].image_x;
				b(2 * index + 1, 0) = aligned[i][j].image_y;
				//printf("(%d %d) (%d %d)\n", aligned2[i][j].image_x, aligned2[i][j].image_y, aligned[i][j].image_x, aligned[i][j].image_y);
				index++;
			}
	Matrix<double> U(6, 6), S(6, 6), VT(6, 6);
	Vector<double> s(6);
	svd(A, U, s, VT);
	for (int i = 0; i < 6; i++)
	{
		if (s(i) > 10e-15) // don't divide by 0
			S(i, i) = 1.0 / s(i);
	}
	Matrix<double> result = transpose(VT) * S * transpose(U) * b;
	Matrix3x3 H;
	H(0, 0) = result(0, 0); H(0, 1) = result(1, 0); H(0, 2) = result(2, 0);
	H(1, 0) = result(3, 0); H(1, 1) = result(4, 0); H(1, 2) = result(5, 0);
	H(2, 0) = 0; H(2, 1) = 0; H(2, 2) = 1;
	return H;
}

int main( int argc, char *argv[] )
{
	// code to composite two images together
	// x = H1^{-1} * trans * H2 y
	/*Matrix3x3 H1(0.999198,0.0020664,41.7295,-0.000307165,1.00242,-47.5016, 0,0,1);
	Matrix3x3 H2(0.999323,-0.00175819,52.1883,-0.000745662,1.00171,-41.3955,0,0,1);
	Matrix3x3 H3(0.995601,-0.00335033,-1229.49,-0.000102755,0.998386,7.31238, 0, 0, 1);
	//Matrix3x3 H3(1, 0, 52.18 - 41.7295 -(2603 - 1364), 0, 1, 5 - 41.39 + 47.5, 0, 0, 1);
	boost::shared_ptr<DiskImageResource> rsrc1(new DiskImageResourceIsis("../data/Apollo15-CAM-CUB/AS15-M-1134.lev1.cub"));
	boost::shared_ptr<DiskImageResource> rsrc2(new DiskImageResourceIsis("../data/Apollo15-CAM-CUB/AS15-M-1135.lev1.cub"));
	DiskImageView<PixelGray<float> > temp1(rsrc1), temp2(rsrc2);

	ImageView<PixelGray<float> > cub1 = normalize(temp1);
	ImageView<PixelGray<float> > cub2 = normalize(temp2);
	cub2 = transform(cub2, HomographyTransform(H3));
	crop(cub1, 0, cub1.cols() / 2 + 100, (cub1.rows() + (int)H3(0, 2)) / 2 + 1100, (cub1.cols() - (int)H3(1, 2)) / 2 - 100) = crop(cub2, 0, cub1.cols() / 2 + 100, (temp1.rows() + (int)H3(0, 2)) / 2 + 1100, (temp1.cols() - (int)H3(1, 2)) / 2 - 100);
	//crop(cub2, cub1.rows() / 2, 0, cub1.rows() / 2, cub1.cols()) = crop(cub1, cub1.rows() / 2, 0, cub1.rows() / 2, cub1.cols());
	ImageView<PixelRGB<uint8> > assembledImg(temp1.rows(), temp1.cols());
	assembledImg = apply_mask(cub1 * 255, 0);
	write_image("mosaic.png", assembledImg);
	return 0;*/


	// command line options
	string inputCSVFilename; 
	
	std::string inputCubFile, tracksListFile, gcpFile;
	std::string outputFile, dataFile, imageFile, startMatricesFilename, pyramidFile;
	bool globalAlignment = false;
	int transSearchWindow = 0, transSearchStep = 0;
	float thetaSearchWindow = 0.0, thetaSearchStep = 0.0;

	po::options_description general_options("Options");
	general_options.add_options()
	("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
	("tracksList,t", po::value<std::string>(&tracksListFile))
	("inputCubFile,i", po::value<std::string>(&inputCubFile))
	("imagePyramid,p", po::value<std::string>(&pyramidFile))
	("outputFile,o", po::value<std::string>(&outputFile))
	("dataFile,d", po::value<std::string>(&dataFile))
	("outputImage", po::value<std::string>(&imageFile))
	("startMatrices,s", po::value<std::string>(&startMatricesFilename))
	("globalAlignment", po::value<bool>(&globalAlignment))
	("transSearchWindow", po::value<int>(&transSearchWindow))
	("transSearchStep", po::value<int>(&transSearchStep))
	("thetaSearchWindow", po::value<float>(&thetaSearchWindow))
	("thetaSearchStep", po::value<float>(&thetaSearchStep))
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
	if (vm.count("startMatrices"))
		trackTransforms = load_track_transforms(startMatricesFilename);
	else
		for (unsigned int i = 0; i < trackPts.size(); i++)
			trackTransforms.push_back(Matrix3x3(1, 0, 0, 0, 1, 0, 0, 0, 1));
	
	// use this later, only if saving to file
	vector<gcp> gcpArray;
	if (gcpFile.length() > 0)
		gcpArray = ComputeSalientLOLAFeatures(trackPts);

	vector<vector< AlignedLOLAShot> > aligned;

	if (vm.count("inputCubFile") > 0)
	{
		aligned = align_to_image(trackPts, inputCubFile, trackTransforms, globalAlignment,
				transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep, imageFile);
	}
	else if (vm.count("imagePyramid") > 0)
	{
		aligned = align_to_image_pyramid(trackPts, pyramidFile, trackTransforms, globalAlignment,
				transSearchWindow, transSearchStep, thetaSearchWindow, thetaSearchStep, imageFile);
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
	
	if (gcpFile.length() > 0)
	{
		UpdateGCP(aligned, "../data/Apollo15-CAM-CUB/AS15-M-1134.lev1.cub", gcpArray);
		SaveGCPoints(gcpArray, gcpFile);
	}

	if (dataFile.length() > 0)
		save_track_data(aligned, dataFile);
 
	return 0;

}

