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
#include "featuresLOLA.h"

void SaveGCPoints(vector<vector<LOLAShot> > trackPts,  std::vector<std::string> DRGFiles,  std::vector<int> overlapIndices, 
                  vector<Vector<float, 6> > optimalTransfArray, vector<float> optimalErrorArray, string gcpFilename)
{


  //for all features in the LOLA data
  //{
  float x = 1;
  float y = 1; 
  float z = 1;
  float sigma_x = 1;
  float sigma_y = 1; 
  float sigma_z = 1;
 
 
  stringstream ss;
  ss<<0;
  cout<<"string_num"<<ss.str();
  string this_gcpFilename = gcpFilename+"_"+ss.str()+".txt";
  cout<<this_gcpFilename<<endl;
  
  FILE *fp = fopen(this_gcpFilename.c_str(), "w");
  fprintf(fp, "%f %f %f %f %f %f\n", x, y, z, sigma_x, sigma_y, sigma_z);
  
  for (int k = 0; k < overlapIndices.size(); k++){
       //int i = (int) floor(finalTransfArray[0]*trackPts[ti][si].imgPt[li].x + finalTransfArray[1]*trackPts[ti][si].imgPt[li].y + finalTransfArray[2]);
       //int j = (int) floor(finalTransfArray[3]*trackPts[ti][si].imgPt[li].x + finalTransfArray[4]*trackPts[ti][si].imgPt[li].y + finalTransfArray[5]);

       float i = optimalTransfArray[k][0]*x + optimalTransfArray[k][1]*y + optimalTransfArray[k][2];
       float j = optimalTransfArray[k][3]*x + optimalTransfArray[k][4]*y + optimalTransfArray[k][5];

       fprintf(fp, "%s %f %f\n", DRGFiles[overlapIndices[k]].c_str(), i, j);
  }
  
  fclose(fp);
  //}
  
}
void printOverlapList(std::vector<int>  overlapIndices)
{
  printf("numOverlapping images = %d\n", (int)(overlapIndices.size()));
    for (int i = 0; i < overlapIndices.size(); i++){
      printf("%d ", overlapIndices[i]);
    }
}

//this will be used to compute the makeOverlapList in a more general way.
//it takes into consideration any set of overlapping images.
Vector4 ComputeGeoBoundary(GeoReference Geo, int width, int height)
{
 
  Vector4 corners;
  Vector2 leftTopPixel(0,0);
  Vector2 leftTopLonLat = Geo.pixel_to_lonlat(leftTopPixel);

  Vector2 rightBottomPixel(width-1, height-1);
  Vector2 rightBottomLonLat = Geo.pixel_to_lonlat(rightBottomPixel);
  
  float minLon = leftTopLonLat(0);
  float minLat = leftTopLonLat(1);
  float maxLon = rightBottomLonLat(0);
  float maxLat = rightBottomLonLat(1);

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

       DiskImageView<PixelMask<PixelGray<uint8> > >  image(inputFiles[i]);
       GeoReference geo;
       read_georeference(geo, inputFiles[i]);
       int width = image.cols();
       int height = image.rows();
       Vector4 corners = ComputeGeoBoundary(geo, width, height);
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

int main( int argc, char *argv[] ) {

  string inputCSVFilename; 
  std::string configFilename="coregister_settings.txt";
  
  std::vector<std::string> DRGFiles;
  std::vector<std::string> DEMFiles;
  std::string resDir = "../results";
 
  po::options_description general_options("Options");
  general_options.add_options()
    ("LidarFile,l", po::value<std::string>(&inputCSVFilename))
    ("DRGFiles,i", po::value<std::vector<std::string> >(&DRGFiles))
    ("res-directory,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
    ("config-filename,c", po::value<std::string>(&configFilename)->default_value("coregister_settings.txt"), "configuration filename.")
    ("help,h", "Display this help message");
  

  po::options_description hidden_options("");

  hidden_options.add_options()
    ("DEMFiles,d", po::value<std::vector<std::string> >(&DEMFiles));
 
  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("DRGFiles", -1);

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

  if(( vm.count("DRGFiles") < 1 ) && (vm.count("DEMFiles") < 1)) {
    std::cerr << "Error: Must specify at least one orthoprojected image file or one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }


  struct CoregistrationParams settings;
  ReadConfigFile((char*)configFilename.c_str(), &settings);
  PrintGlobalParams(&settings);
 
  ModelParams modelParams;
  string modelParamsFilename="light_viewer_pos.txt";
  ReadModelParamsFile(modelParamsFilename, &modelParams);
  PrintModelParams(&modelParams);
  
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
  std::vector<int> overlapIndices = makeOverlapList(DRGFiles, lon_lat_bb);
  printOverlapList(overlapIndices);
  printf("done\n");

  //int numOverlappingImages = (int)overlapIndices.size();
  int numOverlappingImages = 1;
  //TO DO: determine the LOLA features
 
  //allocate memory for all final transforms, one per image
  vector<Vector<float, 6> > optimalTransfArray;
  vector<float> optimalErrorArray;
  optimalTransfArray.resize(numOverlappingImages);
  optimalErrorArray.resize(numOverlappingImages);

  for (int k = 0; k < numOverlappingImages; k++){
    
    string inputDEMFilename;
    string inputDRGFilename = DRGFiles[k];  
    cout<<inputDRGFilename<<endl;

    string DRGFilenameNoPath = sufix_from_filename(inputDRGFilename);

    string imgPtsFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "img.txt";    
    string reflectancePtsFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "refl.txt";  
    string syntImgPtsFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "synth.txt";  
    string altitudePtsFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "alt.txt";
    string demPtsFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "dem.txt";
    string lolaTracksFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "_lola.tif";  
    string outFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "_results.tif";  
    string lolaFeaturesFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "_features_lola.txt";  
    string transformFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "transform_results"  +".txt";  
    string matchResultsFilename = resDir + prefix_less3_from_filename(DRGFilenameNoPath) + "match_results"  +".txt";  

    //create the results directory and prepare the output filenames - END

    vector<Vector<float, 6> >initTransfArray;
    initTransfArray.resize(settings.maxNumStarts);
    for (int i = 0; i < settings.maxNumStarts; i++){
      initTransfArray[i][0] = 1.0;
      initTransfArray[i][1] = 0.0;
      initTransfArray[i][2] = (i-settings.maxNumStarts/2)*5;
      initTransfArray[i][3] = 0.0;
      initTransfArray[i][4] = 1.0;
      initTransfArray[i][5] = 0.0;//(i-maxNumStarts/2)*25;
    }    
 
    vector<int> trackIndices;
    trackIndices.resize(trackPts.size());
    for (int i = 0; i < trackPts.size(); i++){
      trackIndices[i] = i;
    }
  
    vector<Vector<float, 6> > finalTransfArray;
    vector<float> errorArray;
    errorArray.resize(settings.maxNumStarts);
    finalTransfArray.resize(settings.maxNumStarts);

    //initialization step for LIMA - START  
  
    DiskImageView<PixelGray<uint8> >   DRG(inputDRGFilename);
    GeoReference DRGGeo;
    read_georeference(DRGGeo, inputDRGFilename);


    ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
									ConstantEdgeExtension()),
							    BilinearInterpolation());
  
    //get the true image points
    cout << "GetAllPtsFromImage..." << endl; 
    GetAllPtsFromImage(trackPts, interpDRG, DRGGeo);
  
    cout << "ComputeTrackReflectance..." << endl;
    ComputeAllReflectance(trackPts, modelParams, settings);
    //initialization step for LIMA - END  

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

 
    if (settings.useLOLAFeatures){
      cout << "Computing the LOLA features and weights ... ";
      int halfWindow = 10;
      //float topPercent = 0.10;
      ComputeWeights( trackPts, halfWindow, (float)(settings.topPercentFeatures)/100.0, lolaFeaturesFilename);
      cout<<"done."<<endl;
    }

    if (settings.matchingMode == LIMA){
      //return matching error and transform
      cout << "UpdateMatchingParams ..."<< endl;
      UpdateMatchingParamsLIMA_MP(trackPts, inputDRGFilename, 
				  modelParams, settings,  
				  initTransfArray, finalTransfArray, errorArray);
    }

    if (settings.matchingMode == LIDEM){
      //return matching error and transform
      cout << "UpdateMatchingParams ..."<< endl;
      UpdateMatchingParamsLIDEM_MP(trackPts, inputDEMFilename, 
				   modelParams, settings,  
				   initTransfArray, finalTransfArray, errorArray);
    }

    int bestResult = 0;
    float smallestError = 10000000000;
    for (int index = 0; index < initTransfArray.size(); index++){
      printf("OUT %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
	     index, errorArray[index], 
	     finalTransfArray[index](0), finalTransfArray[index](1), 
	     finalTransfArray[index](2), finalTransfArray[index](3),
	     finalTransfArray[index](4), finalTransfArray[index](5));
      if  ( (errorArray[index] < smallestError) && (errorArray[index] > 0) ){
	smallestError = errorArray[index]; 
	bestResult = index;
      }    
    }    
    cout<<"bestResult= "<<bestResult<<endl;
    
    //copy the best transform to optimalTransfArray
    optimalTransfArray[k] = finalTransfArray[bestResult];
    optimalErrorArray[k]  = errorArray[bestResult];

    //write finalTransfArray and errorArray to file
    //this is mostly for debugging info
    SaveReportFile(trackPts, finalTransfArray, errorArray, transformFilename);
    
    //this is what is needed by further steps.
    SaveImagePts(trackPts, finalTransfArray[bestResult], errorArray[bestResult], matchResultsFilename);
  
    //write the image interest point 
    //look at all points s.t. trackPts[si].weight_lsq = 1.0;  
    //determine the image feature location using the finalTransfArray
    //save the image features to file.

    if (settings.displayResults){
      //write results to image outside matching
      ShowFinalTrackPtsOnImage(trackPts, finalTransfArray[bestResult], 
			       trackIndices, inputDRGFilename, outFilename);
    }
  }

 
  //save the GCP
  string gcpFilenameRoot = resDir + prefix_from_filename(sufix_from_filename(inputCSVFilename));
  SaveGCPoints(trackPts, DRGFiles,  overlapIndices, 
               optimalTransfArray, optimalErrorArray, gcpFilenameRoot);


  return 0;
}


















