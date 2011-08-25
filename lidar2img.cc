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


int main( int argc, char *argv[] ) {

  string inputCSVFilename; 
  std::string configFilename="lidar2img_settings.txt";
 
  std::vector<std::string> camCubFiles;
  std::string resDir = "../results";
  std:: string auxDir = "../aux";
  std::string mapCubDir = "../data/map";
  std::string drgDir = "../data/drg";
  std::string camCubFileList;

  po::options_description general_options("Options");
  general_options.add_options()
    ("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
    ("camCubFiles,c", po::value<std::vector<std::string> >(&camCubFiles))
    ("camCubFileList,a", po::value<std::string> (&camCubFileList))
    ("mapCub-directory,m", po::value<std::string>(&mapCubDir)->default_value("../data/map"), "map cub directory.")
    ("drg-directory,d", po::value<std::string>(&drgDir)->default_value("../data/drg"), "drg directory.")
    ("results-directory,r", po::value<std::string>(&resDir)->default_value("../results"), "results directory.")
    ("settings-filename,s", po::value<std::string>(&configFilename)->default_value("lidar2img_settings.txt"), "settings filename.")
    ("help,h", "Display this help message");
  
  po::options_description options("Allowed Options");
  options.add(general_options);

  po::positional_options_description p;
  p.add("camCubFileList", -1);

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

  if(( vm.count("camCubFileList") < 1 )) {
    std::cerr << "Error: Must specify at least one cub image file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }


  //#if 0
  struct CoregistrationParams settings;
  if( ReadConfigFile(configFilename, &settings) ){
      std::cerr << "Config file " << configFilename << " found." << endl;
  }
  else{
    std::cerr << "Config file " << configFilename << " not found, using defaults." << endl;
  }
  //PrintGlobalParams(&settings);
  std::cerr << settings << endl;

 
  ReadCamFileList(camCubFileList, camCubFiles);

  int numCubFiles = camCubFiles.size();
  printf("numCubFiles = %d\n", numCubFiles);

  //create the results directory and prepare the output filenames - START
  system("mkdir ../results"); 
  system("mkdir ../aux");
 
  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);

  //select the overlapping images - START
 
  cout<<"Selecting the overlapping images ..."<<endl;

  std::vector<int> overlapIndices;
  string overlapListFilename = auxDir+string("/")+sufix_from_filename(inputCSVFilename);
  FindAndReplace(overlapListFilename, ".csv", "_overlap_list.txt"); 
  int fileFound = ReadOverlapList(overlapListFilename, overlapIndices);

  if (fileFound == 0){
    Vector4 lat_lon_bb = FindMinMaxLat(trackPts); 
    Vector4 lon_lat_bb;
    lon_lat_bb[0]=lat_lon_bb[2];
    lon_lat_bb[1]=lat_lon_bb[3];
    lon_lat_bb[2]=lat_lon_bb[0];
    lon_lat_bb[3]=lat_lon_bb[1];
    printf("lidar corners: %f %f %f %f\n", lon_lat_bb[0], lon_lat_bb[1], lon_lat_bb[2], lon_lat_bb[3]);
    overlapIndices = makeOverlapList(camCubFiles, lon_lat_bb);
    SaveOverlapList(overlapListFilename, overlapIndices);
  }

  PrintOverlapList(overlapIndices);
  int numOverlappingImages;
  if ((overlapIndices.size() == 1) && (overlapIndices[0] == -1)){
    numOverlappingImages = 0;
  }
  else{
     numOverlappingImages = (int)overlapIndices.size();
  }
  cout<<"numOverlapImages="<<numOverlappingImages<<endl;
  cout<<"done."<<endl;

  //select the overlapping images - END
  
  //ComputeAverageShotDistance(trackPts);
  //ComputeAverageIntraShotDistance(trackPts);
 
  //compute LOLA features - START
  int windowSize = 12;//16
  float salientFeatureThresh = 0.008;//0.015;

  vector<float> filter;
  filter =  MakeLidarFilter(windowSize);
  for (int ti = 0; ti < (int)trackPts.size(); ti++){
    cout<<"trackIndex="<<ti<<" of "<<trackPts.size()<<endl;
    ComputeSalientLOLAFeature(trackPts[ti], filter, salientFeatureThresh);
  }
  //compute LOLA features - END


  cout<<"create the GCPs ..."<<endl;
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
  cout<<"done."<<endl;

 //allocate memory for all final transforms, one per image
  vector<Vector<float, 6> > optimalTransfArray;
  vector<float> optimalErrorArray;
  optimalTransfArray.resize(numOverlappingImages);
  optimalErrorArray.resize(numOverlappingImages);

  cout<<"generating the initTransforms ..."<<endl;
  vector<Vector<float, 6> >initTransfArray;
  GenerateInitTransforms(initTransfArray, settings);
  cout<<"done."<<endl;  

  vector<Vector<float, 6> > bestInitTransfArray; 
  vector<Vector<float, 6> > finalTransfArray;
  vector<float> initMatchingErrorArray;
  vector<float> finalMatchingErrorArray;
  Vector2 transfCentroid;
 
  for (int k = 0; k < numOverlappingImages; k++){
    
    string overlapCamCubFile;
    string overlapMapCubFile;

    overlapCamCubFile = camCubFiles[overlapIndices[k]];
    overlapMapCubFile = GetFilenameNoPath(overlapCamCubFile); 
    overlapMapCubFile = mapCubDir+string("/")+overlapMapCubFile;
    FindAndReplace(overlapMapCubFile, ".cub", "_map.cub"); 

    cout<<"camCubFile="<<overlapCamCubFile<<", mapCubFile="<<overlapMapCubFile<<endl;

    camera::IsisCameraModel model(overlapCamCubFile);
    Vector3 center_of_moon(0,0,0);
    Vector2 pixel_location = model.point_to_pixel( center_of_moon );
    Vector3 cameraPosition = model.camera_center( pixel_location );
    Vector3 lightPosition = model.sun_position( pixel_location );

    //initialization step for LIMA - START  
    cout<<"GetAllPtsFromCub "<<endl; 
    int numValidImgPts = GetAllPtsFromCub(trackPts, overlapMapCubFile);
    cout<<"numValidImgPts="<<numValidImgPts<<endl;
    cout<<"done."<<endl; 
  
    cout << "ComputeTrackReflectance..." << endl;
    int  numValidReflPts = ComputeAllReflectance(trackPts,  cameraPosition, lightPosition, settings);
    cout<<"numValidReflPts="<<numValidReflPts<<endl;
    cout<<"done."<<endl;
    //initialization step for LIMA - END 

  
    if (numValidReflPts > 100){
      //param init LIMA - START 


      vector<Vector4> matchArray = FindMatches2D(trackPts, overlapMapCubFile, settings.matchWindowHalfSize, 40, initMatchingErrorArray);
      
      /*  
      EstimateAffineTransform(bestInitTransfArray, trackPts, matchArray, initMatchingErrorArray);
     
      //EstimateMatchingParamsFromCub(trackPts, overlapMapCubFile, settings.matchWindowHalfSize, bestInitTransfArray, initMatchingErrorArray);
      
      // cout<<"Initializing the affine tranformation ..."<<endl;
      //InitMatchingParamsFromCub(trackPts, overlapMapCubFile,  settings,  
      //			  initTransfArray, bestInitTransfArray, initMatchingErrorArray);
    

      //param init LIMA - END
      
      //iterative matching step for LIMA - START
      int refinedMatching = 0;//1;
      
      if (refinedMatching == 1){
	cout<<"Computing the affine tranformation ..."<<endl;
        
	UpdateMatchingParamsFromCub(trackPts, overlapMapCubFile,  settings.maxNumIter,  
				    bestInitTransfArray, initMatchingErrorArray,
				    finalTransfArray, finalMatchingErrorArray, transfCentroid);
	
        cout<<"done."<<endl;
      }
      else{
	finalTransfArray = bestInitTransfArray;
	finalMatchingErrorArray = initMatchingErrorArray;
      }
      
      GetBestTransform(finalTransfArray, finalMatchingErrorArray, optimalTransfArray[k], optimalErrorArray[k]);
      //iterative matching step for LIMA - END
      */
      
      //save to GCP structure. 
      //UpdateGCP(trackPts, optimalTransfArray[k], overlapCamCubFile, overlapMapCubFile, gcpArray, transfCentroid, 4.0);
      UpdateGCP(trackPts, matchArray, initMatchingErrorArray, overlapCamCubFile, overlapMapCubFile, gcpArray, 4.0);
      
      matchArray.clear();
    }
  }  
  //end matching
 
  //save the GCP
  //string gcpFilenameRoot = resDir + prefix_from_filename(sufix_from_filename(inputCSVFilename));  
  string gcpFilenameRoot = resDir+"/"+GetFilenameNoExt(GetFilenameNoPath(inputCSVFilename));


  cout<<"writting the GCP file..."<<endl;
  SaveGCPoints(gcpArray,  gcpFilenameRoot);
  cout<<"done"<<endl;

  //these are the display functions
  if (settings.analyseFlag){
 
    for (int k = 0; k < numOverlappingImages; k++){
      
      string overlapCamCubFile = camCubFiles[overlapIndices[k]];
      string overlapCamCubFileNoPath = GetFilenameNoPath(overlapCamCubFile); 
 
      string overlapMapCubFile = mapCubDir+string("/")+overlapCamCubFileNoPath;
      FindAndReplace(overlapMapCubFile, ".cub", "_map.cub"); 
     
      string overlapDRGFilename = drgDir+string("/")+overlapCamCubFileNoPath;
      FindAndReplace(overlapDRGFilename, ".cub", "_drg.tif");
 
      string thisResDir = resDir+string("/")+overlapCamCubFileNoPath;
      FindAndReplace(thisResDir, ".cub", "");
      string command = string("mkdir ")+thisResDir;
      system(command.c_str()); 

      //camera::IsisCameraModel model(overlapCamCubFile);
      camera::IsisCameraModel model(overlapMapCubFile);
      Vector3 center_of_moon(0,0,0);
      Vector2 pixel_location = model.point_to_pixel( center_of_moon );
      Vector3 cameraPosition = model.camera_center( pixel_location );
      Vector3 lightPosition = model.sun_position( pixel_location );
     
      string overlapCamCubFileNoPathNoExt = GetFilenameNoExt(overlapCamCubFileNoPath);
      string LOLAFilenameNoPathNoExt = GetFilenameNoExt(GetFilenameNoPath(inputCSVFilename));
  
      string altitudePtsFilename = thisResDir + "/" + LOLAFilenameNoPathNoExt + "_alt.txt";  
      string imgPtsFilename = thisResDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapCamCubFileNoPathNoExt + "_img.txt";  
      string drgPtsFilename = thisResDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapCamCubFileNoPathNoExt + "_drg.txt";  
      string reflectancePtsFilename = thisResDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapCamCubFileNoPathNoExt + "_refl.txt";  
      string syntImgPtsFilename = thisResDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapCamCubFileNoPathNoExt + "_synth.txt";  
       
      //initialization step for LIMA - START  
      //cout<<"GetAllPtsFromCub "<<endl; 
      //int numValidImgPts = GetAllPtsFromCub(trackPts, overlapMapCubFile);
      //cout<<"numValidImgPts="<<numValidImgPts<<endl;
      //cout<<"done."<<endl; 

      cout<<"GetAllPtsFromImage "<<endl; 
      cout<<"overlapDRGFilename="<<overlapDRGFilename<<endl;
      boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(overlapDRGFilename) );
      DiskImageView<PixelGray<uint8> > DRG( rsrc );
      GeoReference DRGGeo;
      read_georeference(DRGGeo, overlapDRGFilename);
      GetAllPtsFromImage(trackPts, DRG, DRGGeo);
      cout<<"done."<<endl; 

 
  
      cout << "ComputeTrackReflectance..." << endl;
      cout<<cameraPosition<<endl;
      cout<<lightPosition<<endl;
      int  numValidReflPts = ComputeAllReflectance(trackPts,  cameraPosition, lightPosition, settings);
      cout<<"numValidReflPts="<<numValidReflPts<<endl;
      cout<<"done."<<endl;
      //initialization step for LIMA - END 
      
      Vector2 gain_bias = ComputeGainBiasFactor(trackPts);
      //SaveImagePoints(trackPts, 1, imgPtsFilename);
      SaveImagePoints(trackPts, 1, drgPtsFilename);
      SaveAltitudePoints(trackPts, 1, altitudePtsFilename);
      SaveReflectancePoints(trackPts, 1.0, reflectancePtsFilename);
      SaveReflectancePoints(trackPts, gain_bias, syntImgPtsFilename);
     
      /*
      string lolaInitTracksFilename = resDir + prefix_less3_from_filename(overlapCamCubFileNoPath) + "init_lola.tif";  
      string lolaInitTracksOnImageFilename = resDir + prefix_less3_from_filename(overlapCamCubFileNoPath) + "init_img_lola.tif";  
      string lolaFinalTracksOnImageFilename = resDir + prefix_less3_from_filename(overlapCamCubFileNoPath) + "final_img_lola.tif"; 
      vector<int> trackIndices;
      trackIndices.resize(trackPts.size());
      for (unsigned int i = 0; i < trackPts.size(); i++){
	   trackIndices[i] = i;
      }

      Vector<float, 6> unitTransfArray;
      unitTransfArray[0]=1;
      unitTransfArray[1]=0;
      unitTransfArray[2]=0;
      unitTransfArray[3]=0;
      unitTransfArray[4]=1;
      unitTransfArray[5]=0;
      ShowFinalTrackPtsOnImage(trackPts, unitTransfArray, 
			       trackIndices, overlapMapCubFile, lolaInitTracksOnImageFilename);
      

      //write results to image outside matching
      //ShowFinalTrackPtsOnImage(trackPts, optimalTransfArray[k], 
      //			 trackIndices, mapCubFiles[k], lolaFinalTracksOnImageFilename);
      */
    }

    //TO DO: this should become all one function - START  
    cout<<"saving the GCP images..."<<endl;
    int gc_index = 0;
 
    for (unsigned int t=0; t<trackPts.size(); t++){
      int featureIndex = 0;
      for (unsigned int s=0; s<trackPts[t].size(); s++){
	if (trackPts[t][s].featurePtLOLA==1){ 
	  if (gcpArray[gc_index].filename.size() > 0){ //non-empty GCP
	    
            stringstream trackIndexString;
            trackIndexString<<gcpArray[gc_index].trackIndex;            

            stringstream featureIndexString;
            featureIndexString<<gcpArray[gc_index].shotIndex;

	    string gcpFilename = gcpFilenameRoot+"_"+trackIndexString.str()+"_"+featureIndexString.str()+".gcp";
	    string assembledImgFilename = gcpFilenameRoot+"_img_"+trackIndexString.str()+"_"+featureIndexString.str()+".tif";
	    cout<<"gcpFilename="<<gcpFilename<<endl;
	    cout<<"assembledImgFilename="<<assembledImgFilename<<endl;
	    SaveGCPImages(gcpArray[gc_index], assembledImgFilename);
	    
	  }
	  featureIndex++;
          gc_index++;
	}
      }
    }

    cout<<"done."<<endl;
    //TO DO: this should become all one function - END  

  }
  //#endif
 
  return 0;

}


















