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
#include <sstream>

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
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;
using namespace vw::camera;
using namespace std;

#include <math.h>
#include "tracks.h"
#include "io.h"
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"
#include "featuresLOLA.h"

void printOverlapList(std::vector<int>  overlapIndices)
{
  printf("numOverlapping images = %d\n", (int)(overlapIndices.size()));
    for (int i = 0; i < overlapIndices.size(); i++){
      printf("%d ", overlapIndices[i]);
    }
}


//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
Vector4 ComputeGeoBoundary(string cubFilename)
{
  GeoReference moonref( Datum("D_MOON"), identity_matrix<3>() );
  boost::shared_ptr<IsisCameraModel> isiscam( new IsisCameraModel( cubFilename ) );
  BBox2 camera_boundary =
    camera_bbox( moonref, boost::shared_dynamic_cast<CameraModel>(isiscam),
		 isiscam->samples(), isiscam->lines() );
 
  Vector4 corners;
 

  float minLon = camera_boundary.min()[0];
  float minLat = camera_boundary.min()[1];
  float maxLon = camera_boundary.max()[0];
  float maxLat = camera_boundary.max()[1];

  printf("minLon = %f, minLat = %f, maxLon = %f, maxLat = %f\n", minLon, minLat, maxLon, maxLat);
  if (maxLat<minLat){
      float temp = minLat;
      minLat = maxLat;
      maxLat = temp;    
  }

  if (maxLon<minLon){
     float temp = minLon;
     minLon = maxLon;
     maxLon = temp;    
  }

  corners(0) = minLon;
  corners(1) = maxLon;
  corners(2) = minLat;
  corners(3) = maxLat;

  return corners;
}

//this function determines the image overlap for the general case
//it takes into consideration any set of overlapping images.
std::vector<int> makeOverlapList(std::vector<std::string> inputFiles, Vector4 currCorners) {
  
  std::vector<int> overlapIndices;
 
  for (unsigned int i = 0; i < inputFiles.size(); i++){

       int lonOverlap = 0;
       int latOverlap = 0; 
     
       Vector4 corners = ComputeGeoBoundary(inputFiles[i]);

       printf("corners = %f %f %f %f\n", corners[0], corners[1], corners[2], corners[3]);       

       if(  ((corners(0)>currCorners(0)) && (corners(0)<currCorners(1))) //minlon in interval 
	    ||((corners(1)>currCorners(0)) && (corners(1)<currCorners(1)))) //maxlon in interval
       {
         lonOverlap = 1;
       }
       if(  ((corners(2)>currCorners(2)) && (corners(2)<currCorners(3))) //minlat in interval 
	    ||((corners(3)>currCorners(2)) && (corners(3)<currCorners(3)))) //maxlat in interval
       {
         latOverlap = 1;
       }
     
       if ((lonOverlap == 1) && (latOverlap == 1)){
           overlapIndices.push_back(i);
       }
  }

  return overlapIndices;
}

void GetBestTransform(vector<Vector<float, 6> > &finalTransfArray,  vector<float> &finalMatchingErrorArray, 
                      Vector<float, 6> &optimalTransf, float &optimalError)
{
    //this will be made a special function to make sure it is flexible wrt to the name of the parameters.
    int bestResult = 0;
    float smallestError = std::numeric_limits<float>::max();
    for (int index = 0; index < finalTransfArray.size(); index++){
      printf("refined %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
	     index, finalMatchingErrorArray[index], 
	     finalTransfArray[index](0), finalTransfArray[index](1), 
	     finalTransfArray[index](2), finalTransfArray[index](3),
	     finalTransfArray[index](4), finalTransfArray[index](5));
      if  ( (finalMatchingErrorArray[index] < smallestError) && (finalMatchingErrorArray[index] > 0) ){
	smallestError = finalMatchingErrorArray[index]; 
	bestResult = index;
      }    
    } 
    cout<<"bestResult= "<<bestResult<<endl;

    //copy the best transform to optimalTransfArray
    optimalTransf = finalTransfArray[bestResult];
    optimalError  = finalMatchingErrorArray[bestResult];
}

void GenerateInitTransforms( vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings)
{
    initTransfArray.resize(settings.maxNumStarts);
    for (int i = 0; i < settings.maxNumStarts; i++){
      initTransfArray[i][0] = 1.0;
      initTransfArray[i][1] = 0.0;
      initTransfArray[i][2] = (i-settings.maxNumStarts/2)*5;
      initTransfArray[i][3] = 0.0;
      initTransfArray[i][4] = 1.0;
      initTransfArray[i][5] = 0.0;//(i-maxNumStarts/2)*25;
    }  
}

int main( int argc, char *argv[] ) {

  string inputCSVFilename; 
  std::string configFilename="coregister_settings.txt";
  
  //std::vector<std::string> DRGFiles;
  std::vector<std::string> DEMFiles;
  std::vector<std::string> cubFiles;
  std::string resDir = "../results";
 
  po::options_description general_options("Options");
  general_options.add_options()
    ("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
    ("cubFiles,c", po::value<std::vector<std::string> >(&cubFiles))
    ("results-directory,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
    ("settings-filename,s", po::value<std::string>(&configFilename)->default_value("coregister_settings.txt"), "settings filename.")
    ("help,h", "Display this help message");
  

  po::options_description hidden_options("");

  hidden_options.add_options()
    ("DEMFiles,d", po::value<std::vector<std::string> >(&DEMFiles));
 
  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("cubFiles", -1);

  std::ostringstream usage;
  usage << "Description: main code for Lidar to image or DEM co-registration" << std::endl << std::endl;
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

  if(( vm.count("cubFiles") < 1 ) /*&& (vm.count("DEMFiles") < 1)*/) {
    std::cerr << "Error: Must specify at least one orthoprojected image file or one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }


  struct CoregistrationParams settings;
  ReadConfigFile((char*)configFilename.c_str(), &settings);
  PrintGlobalParams(&settings);


  int numCubFiles = cubFiles.size();
  printf("numCubFiles = %d\n", numCubFiles);
  vector<ModelParams> modelParamsArray;
  modelParamsArray.resize(numCubFiles);
  for (int i = 0; i < numCubFiles; i++){
    printf("cubFiles[%d] = %s\n", i, cubFiles[i].c_str());
    camera::IsisCameraModel model(cubFiles[i]);
    Vector3 center_of_moon(0,0,0);
    Vector2 pixel_location = model.point_to_pixel( center_of_moon );
    Vector3 cameraPosition = model.camera_center( pixel_location );
    Vector3 lightPosition= model.sun_position( pixel_location );
    modelParamsArray[i].sunPosition[0] = lightPosition[0];
    modelParamsArray[i].sunPosition[1] = lightPosition[1];
    modelParamsArray[i].sunPosition[2] = lightPosition[2];

    modelParamsArray[i].spacecraftPosition[0] = cameraPosition[0];
    modelParamsArray[i].spacecraftPosition[1] = cameraPosition[1];
    modelParamsArray[i].spacecraftPosition[2] = cameraPosition[2];
    
    PrintModelParams(&modelParamsArray[i]);
  }

  //create the results directory and prepare the output filenames - START
  system("mkdir ../results");
  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);
 
  //select the overlapping images
  printf("Selecting the overlapping images ...\n");
  Vector4 lat_lon_bb = FindMinMaxLat(trackPts); 
  Vector4 lon_lat_bb;
  lon_lat_bb[0]=lat_lon_bb[2];
  lon_lat_bb[1]=lat_lon_bb[3];
  lon_lat_bb[2]=lat_lon_bb[0];
  lon_lat_bb[3]=lat_lon_bb[1];
  printf("lidar corners: %f %f %f %f\n", lon_lat_bb[0], lon_lat_bb[1], lon_lat_bb[2], lon_lat_bb[3]);
  
  std::vector<int> overlapIndices = makeOverlapList(cubFiles, lon_lat_bb);
  printOverlapList(overlapIndices);
  printf("done\n");

  int numOverlappingImages = (int)overlapIndices.size();
  //determine the LOLA features
  int halfWindow = 10;
  printf("numTracks = %d\n", (int)trackPts.size());
  for (int ti = 0; ti < trackPts.size(); ti++){
    printf("trackIndex = %d\n", ti);
    ComputeSalientLOLAFeature(trackPts[ti], halfWindow, (float)(settings.topPercentFeatures)/1000.0);
  }

 //allocate memory for all final transforms, one per image
  vector<Vector<float, 6> > optimalTransfArray;
  vector<float> optimalErrorArray;
  optimalTransfArray.resize(numOverlappingImages);
  optimalErrorArray.resize(numOverlappingImages);

  for (int k = 0; k < numOverlappingImages; k++){

    //string inputDEMFilename;
    string inputImgFilename = cubFiles[k];;  
    cout<<inputImgFilename<<endl;

    string imgFilenameNoPath = sufix_from_filename(inputImgFilename);

    string imgPtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "img.txt";    
    string reflectancePtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "refl.txt";  
    string syntImgPtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "synth.txt";  
    string altitudePtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "alt.txt";
    string demPtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "dem.txt";
    string lolaTracksFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "_lola.tif";  
    string outFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "_results.tif";  
    string lolaFeaturesFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "_features_lola.txt";  
    string transformFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "transform_results"  +".txt";  
    string matchResultsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "match_results"  +".txt";  

    //create the results directory and prepare the output filenames - END

    vector<Vector<float, 6> >initTransfArray;
    vector<Vector<float, 6> >bestInitTransfArray; 
    vector<Vector<float, 6> > finalTransfArray;
    vector<float> initMatchingErrorArray;
    vector<float> finalMatchingErrorArray;

    GenerateInitTransforms(initTransfArray, settings);
    
    vector<int> trackIndices;
    trackIndices.resize(trackPts.size());
    for (int i = 0; i < trackPts.size(); i++){
      trackIndices[i] = i;
    }

    //initialization step for LIMA - START  
    cout<<"GetAllPtsFromCub"<<endl; 
    GetAllPtsFromCub(trackPts, cubFiles[k]);

    cout << "ComputeTrackReflectance..." << endl;
    ComputeAllReflectance(trackPts, modelParamsArray[k], settings);
    //initialization step for LIMA - END  

    /*
    if (settings.analyseFlag == 1){
      float scaleFactor = ComputeScaleFactor(trackPts);
      SaveImagePoints(trackPts, 3, imgPtsFilename);
      SaveAltitudePoints(trackPts, 3, altitudePtsFilename);
      SaveReflectancePoints(trackPts, 1.0, reflectancePtsFilename);
      SaveReflectancePoints(trackPts, scaleFactor, syntImgPtsFilename);
      SaveDEMPoints(trackPts, inputDEMFilename, demPtsFilename);

      int numVerPts = 6000;
      int numHorPts = 6000;
      MakeGrid(trackPts, numVerPts, numHorPts, lolaTracksFilename, trackIndices);
    }
    */
 
    if (settings.useLOLAFeatures){
      cout << "Computing the LOLA features and weights ... ";
      int halfWindow = 10; //this should go into settings
      ComputeWeights( trackPts, halfWindow, (float)(settings.topPercentFeatures)/1000.0, lolaFeaturesFilename);
      cout<<"done."<<endl;
    }
   
    cout<<"Initializing the affine tranformation ..."<<endl;
    //initMatchingErrorArray retains the matching error for each bestInitTransfArray
    InitMatchingParamsFromCub(trackPts, cubFiles[k], modelParamsArray[k], settings,  
			      initTransfArray, bestInitTransfArray, initMatchingErrorArray);
    cout<<"done."<<endl;

    cout<<"Computing the affine tranformation ..."<<endl;
    //finalMatchingErrorArray retains the matching error after the last iteration for each finalTransfArray   
    UpdateMatchingParamsFromCub(trackPts, cubFiles[k], modelParamsArray[k], settings.maxNumIter,  
      			        bestInitTransfArray, finalTransfArray, finalMatchingErrorArray);
    cout<<"done."<<endl;

    GetBestTransform(finalTransfArray, finalMatchingErrorArray, optimalTransfArray[k], optimalErrorArray[k]);
    
    if (settings.displayResults){
       //write results to image outside matching
       ShowFinalTrackPtsOnImage(trackPts, optimalTransfArray[k], 
			       trackIndices, cubFiles[k], outFilename);
    }    

  }
 
  //save the GCP
  string gcpFilenameRoot = resDir + prefix_from_filename(sufix_from_filename(inputCSVFilename));

  for (int k = 0; k < numOverlappingImages; k++){
    printf("k=%d, [0]:%f  [1]:%f [2]:%f [3]:%f [4]:%f [5]:%f\n",
	   k, optimalTransfArray[k][0], optimalTransfArray[k][1],
              optimalTransfArray[k][2], optimalTransfArray[k][3],
              optimalTransfArray[k][4], optimalTransfArray[k][5]);
  }
  cout<<"writting the GC file..."<<endl;
  SaveGCPoints(trackPts, cubFiles,  overlapIndices, 
               optimalTransfArray, optimalErrorArray, gcpFilenameRoot);
  return 0;

}


















