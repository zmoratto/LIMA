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
#include "util.h"
#include "tracks.h"
#include "match.h"
#include "coregister.h"
#include "display.h"
#include "weights.h"
#include "featuresLOLA.h"


int main( int argc, char *argv[] ) {

  string inputCSVFilename; 
  std::string configFilename="lidar2img_settings.txt";
  
  std::vector<std::string> DEMFiles;
  std::vector<std::string> mapCubFiles;
  std::string resDir = "../results";
 
  po::options_description general_options("Options");
  general_options.add_options()
    ("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
    ("mapCubFiles,c", po::value<std::vector<std::string> >(&mapCubFiles))
    ("results-directory,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
    ("settings-filename,s", po::value<std::string>(&configFilename)->default_value("lidar2img_settings.txt"), "settings filename.")
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

  if(( vm.count("cubFiles") < 1 )) {
    std::cerr << "Error: Must specify at least one orthoprojected image file or one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }

  //#if 0
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


  int numCubFiles = mapCubFiles.size();
  printf("numCubFiles = %d\n", numCubFiles);
  vector<ModelParams> modelParamsArray;
  modelParamsArray.resize(numCubFiles);
  for (int i = 0; i < numCubFiles; i++){
    printf("mapCubFiles[%d] = %s\n", i, mapCubFiles[i].c_str());
    camera::IsisCameraModel model(mapCubFiles[i]);
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
  
  std::vector<int> overlapIndices = makeOverlapList(mapCubFiles, lon_lat_bb);
  printOverlapList(overlapIndices);
  printf("done\n");

  int numOverlappingImages = (int)overlapIndices.size();
  ComputeAverageShotDistance(trackPts);
  ComputeAverageIntraShotDistance(trackPts);
  /*
  //determine the LOLA features
  int halfWindow = 10;
  printf("numTracks = %d\n", (int)trackPts.size());
  for (int ti = 0; ti < (int)trackPts.size(); ti++){
    printf("trackIndex = %d\n", ti);
    ComputeSalientLOLAFeature(trackPts[ti], halfWindow, (float)(settings.topPercentFeatures)/1000.0);
    ComputeSalientLOLAFeature(vector<LOLAShot > & trackPts, vector<float> filter, float salientFeatureThresh);
  }
  */
  //compute LOLA features
  printf("numTracks = %d\n", (int)trackPts.size());
  /* 
  //build step filter - START
  vector<float> filter;

  int windowSize = 21;
  int halfWindowSize = windowSize/2;
  filter.resize(windowSize);
  for (int i = 0; i < windowSize ;i++ ){
    if (i < halfWindowSize){
      filter[i] = -1.0; 
    }
    if(i == halfWindowSize){
      filter[i] = 0.0;
    }
    if (i > halfWindowSize){
      filter[i] = 1.0;
    }
  }
  //build step filter - END
  float salientFeatureThresh = 0.18;
  */

  //build crater filter - START
  vector<float> filter;
  int windowSize = 12;//16
  int halfWindowSize = windowSize/2;
  int quartWindowSize = windowSize/4;
  filter.resize(windowSize);
  for (int i = 0; i < windowSize ;i++ ){
    if ((i < quartWindowSize) || (i > windowSize-quartWindowSize-1)){
      filter[i] = 1.0; 
    }
    else{
      filter[i] = -1.0; 
    }
    //printf("%d %f\n", i, filter[i]);
  }
  
  //build crater filter -  END
  float salientFeatureThresh = 0.008;//0.015;

  for (int ti = 0; ti < (int)trackPts.size(); ti++){
    printf("trackIndex = %d\n", ti);
    ComputeSalientLOLAFeature(trackPts[ti], filter, salientFeatureThresh);
  }

  //this should be returned by ComputeSalientLOLAFeature
  vector<gcp> gcpArray;
  for (unsigned int t=0; t<trackPts.size(); t++){
    for (unsigned int s=0; s<trackPts[t].size(); s++){
      if (trackPts[t][s].featurePtLOLA==1){
        gcp this_gcp;
	this_gcp.lon = trackPts[t][s].LOLAPt[2].x();
	this_gcp.lat = trackPts[t][s].LOLAPt[2].y(); 
	this_gcp.rad = trackPts[t][s].LOLAPt[2].z()*1000;
        this_gcp.sigma_lon = 1.0;
        this_gcp.sigma_lat = 1.0;
        this_gcp.sigma_rad = 1.0;
        gcpArray.push_back(this_gcp);
      }
    }
  }

 //allocate memory for all final transforms, one per image
  vector<Vector<float, 6> > optimalTransfArray;
  vector<float> optimalErrorArray;
  optimalTransfArray.resize(numOverlappingImages);
  optimalErrorArray.resize(numOverlappingImages);


  for (int k = 0; k < numOverlappingImages; k++){

    //string inputDEMFilename;
    string inputImgFilename = mapCubFiles[k];  
    cout<<inputImgFilename<<endl;

    string imgFilenameNoPath = sufix_from_filename(inputImgFilename);

    string imgPtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "img.txt";    
    string reflectancePtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "refl.txt";  
    string syntImgPtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "synth.txt";  
    string altitudePtsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "alt.txt";
    string lolaFeaturesFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "_features_lola.txt";  
    string transformFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "transform_results"  +".txt";  
    string matchResultsFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "match_results"  +".txt";  
    string lolaTracksFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "lola.tif";  
    string lolaInitTracksOnImageFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "img_lola.tif";  
    string lolaFinalTracksOnImageFilename = resDir + prefix_less3_from_filename(imgFilenameNoPath) + "results_img_lola.tif";  
  
    //create the results directory and prepare the output filenames - END

    vector<int> trackIndices;
    trackIndices.resize(trackPts.size());
    for (int i = 0; i < trackPts.size(); i++){
      trackIndices[i] = i;
    }
    
    if (settings.analyseFlag == 1){
      
      //float scaleFactor = ComputeScaleFactor(trackPts);
      //SaveImagePoints(trackPts, 3, imgPtsFilename);
      //SaveAltitudePoints(trackPts, 3, altitudePtsFilename);
      
      //SaveReflectancePoints(trackPts, 1.0, reflectancePtsFilename);
      //SaveReflectancePoints(trackPts, scaleFactor, syntImgPtsFilename);
      //int numVerPts = 6000;
      //int numHorPts = 6000;
      //MakeGrid(trackPts, numVerPts, numHorPts, lolaTracksFilename, trackIndices); 
      
    }
   
    //initialization step for LIMA - START  
    cout<<"GetAllPtsFromCub"<<endl; 
    GetAllPtsFromCub(trackPts, mapCubFiles[k]);
    
    
    if (settings.analyseFlag == 1){
        Vector<float, 6> unitTransfArray;
        unitTransfArray[0]=1;
        unitTransfArray[1]=0;
	unitTransfArray[2]=0;
	unitTransfArray[3]=0;
	unitTransfArray[4]=1;
	unitTransfArray[5]=0;
        ShowFinalTrackPtsOnImage(trackPts, unitTransfArray, 
			         trackIndices, mapCubFiles[k], lolaInitTracksOnImageFilename);
    }
    

    cout << "ComputeTrackReflectance..." << endl;
    ComputeAllReflectance(trackPts, modelParamsArray[k], settings);
    //initialization step for LIMA - END  

     
    if (settings.useReflectanceFeatures){
      cout << "Computing LOLA reflectance features and weights ... ";
      int halfWindow = 10; //this should go into settings
      ComputeWeights( trackPts, halfWindow, (float)(settings.topPercentFeatures)/1000.0, lolaFeaturesFilename);
      cout<<"done."<<endl;
    }
    else{
      ResetWeights(trackPts);
    }
   

    //start matching
    vector<Vector<float, 6> >initTransfArray;
    vector<Vector<float, 6> >bestInitTransfArray; 
    vector<Vector<float, 6> > finalTransfArray;
    vector<float> initMatchingErrorArray;
    vector<float> finalMatchingErrorArray;
    Vector2 transfCentroid;
    
    GenerateInitTransforms(initTransfArray, settings);

    cout<<"Initializing the affine tranformation ..."<<endl;
    //initMatchingErrorArray retains the matching error for each bestInitTransfArray
    InitMatchingParamsFromCub(trackPts, mapCubFiles[k], modelParamsArray[k], settings,  
			      initTransfArray, bestInitTransfArray, initMatchingErrorArray);
    cout<<"done."<<endl;
   
    int refinedMatching = 1;

    if (refinedMatching == 1){
      cout<<"Computing the affine tranformation ..."<<endl;
      //finalMatchingErrorArray retains the matching error after the last iteration for each finalTransfArray   
      UpdateMatchingParamsFromCub(trackPts, mapCubFiles[k], modelParamsArray[k], settings.maxNumIter,  
				  bestInitTransfArray, finalTransfArray, finalMatchingErrorArray, transfCentroid);
      cout<<"done."<<endl;
    }
    else{
       finalTransfArray = bestInitTransfArray;
       finalMatchingErrorArray = initMatchingErrorArray;
    }

    GetBestTransform(finalTransfArray, finalMatchingErrorArray, optimalTransfArray[k], optimalErrorArray[k]);

    //save to GCP structure.
    cout<<"CENTROID"<<transfCentroid<<endl;
    UpdateGCP(trackPts, optimalTransfArray[k], mapCubFiles[k], gcpArray, transfCentroid);

    if (settings.analyseFlag){
       //write results to image outside matching
       ShowFinalTrackPtsOnImage(trackPts, optimalTransfArray[k], 
			        trackIndices, mapCubFiles[k], lolaFinalTracksOnImageFilename);
    }    
    
  }  
  //end matching
 
  //save the GCP
  string gcpFilenameRoot = resDir + prefix_from_filename(sufix_from_filename(inputCSVFilename));
  /*
  for (int k = 0; k < numOverlappingImages; k++){
    printf("k=%d, [0]:%f  [1]:%f [2]:%f [3]:%f [4]:%f [5]:%f\n",
	   k, optimalTransfArray[k][0], optimalTransfArray[k][1],
              optimalTransfArray[k][2], optimalTransfArray[k][3],
              optimalTransfArray[k][4], optimalTransfArray[k][5]);
  }
  */
  cout<<"writting the GC file..."<<endl;
  SaveGCPoints(gcpArray,  gcpFilenameRoot);


  if (settings.analyseFlag){
 
    //TO DO: this should become all one function - START  
    cout<<"writting the GCP images..."<<endl;
    int gc_index = 0;
  
    for (unsigned int t=0; t<trackPts.size(); t++){
      for (unsigned int s=0; s<trackPts[t].size(); s++){
	if (trackPts[t][s].featurePtLOLA==1){ 
	  if (/*(25*(gc_index/25)==gc_index) &&*/ (gcpArray[gc_index].filename.size() > 0)){ //non-empty GCP
	    stringstream ss;
	    ss<<gc_index;
	    string gcpFilename = gcpFilenameRoot+"_"+ss.str()+".gcp";
	    string assembledImgFilename = gcpFilenameRoot+"_img_"+ss.str()+".tif";
	    cout<<"gcpFilename="<<gcpFilename<<endl;
	    cout<<"assembledImgFilename="<<assembledImgFilename<<endl;
	    SaveGCPImages(gcpArray[gc_index], assembledImgFilename);
	    
	  }
	  gc_index++;
	}
      }
    }
    //TO DO: this should become all one function - END  

  }
  //#endif
  /*
  string gcpFilename = string("../results/RDR_1E8E_21N28NPointPerRow_csv_table_0.gcp");
  string assembledImgFilename = string("../results/RDR_1E8E_21N28NPointPerRow_csv_table_0.tif");
  string cubDirname = string("../data/Apollo15-CUB");
  SaveGCPImages(gcpFilename, cubDirname, assembledImgFilename);
  */
  return 0;

}


















