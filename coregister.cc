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
#include "io.h"
#include "match.h"
#include "coregister.h"
#include "display.h"



void WhenBuildImgPts(vector<LOLAShot> trackPts){
  //this function is to understand which pts GetTrackPtsFromImage will extract

  /* First set of tests:
     1. what is the length of trackPts[0] // [0] is the index being fed in
     2. what is the # accepted based on the id == 3 condition
     3. what is the # accepted on the valid refl. condition?
   */

  //1. # of track points
  //printf("\n\n# trackPts = %d\n",trackPts.size());

  //2. # accepted on id==3
  int number_pass = 0;
  int ID = 3;
  for(int i = 0; i < trackPts.size(); i++){
    for(int k = 0; k < trackPts[i].LOLAPt.size(); k++){

      int id = trackPts[i].LOLAPt[k].s; 
      if (id == ID){
        //  imgPts.push_back((float)DRGVal);
        number_pass += 1;
      }

    }
  }
  printf("number img pts = %d\n",number_pass);

  //3. # valid refl.
  int numb_refl_valid = 0;
  vector<pointCloud> ptHere;
  pointCloud centerPt;
  pointCloud topPt;
  pointCloud leftPt;

  for(int i = 1; i < trackPts.size(); i++){
    ptHere = trackPts[i].LOLAPt;
    centerPt = GetPointFromIndex(ptHere, 3);
    topPt = GetPointFromIndex(ptHere, 2);
    leftPt = GetPointFromIndex(ptHere, 1);

    if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1))
    {
      numb_refl_valid += 1;
    }
  }
  printf("number valid refl pts = %d \n\n",numb_refl_valid);
}


vector<float> AllTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename)
{ 
  DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  vector<float> allImgPts;
  allImgPts.resize(trackPts.size());  
  vector<pointCloud> ptHere;

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  for(int i = 0; i < trackPts.size(); i++){
    ptHere = trackPts[i].LOLAPt;
    pointCloud centerPt  = GetPointFromIndex( ptHere, 3);
    pointCloud topPt     = GetPointFromIndex( ptHere, 2);
    pointCloud leftPt    = GetPointFromIndex( ptHere, 1);

    if((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){

      float lon = centerPt.coords[0];
      float lat = centerPt.coords[1];
      float rad = centerPt.coords[2];

      Vector2 DEM_lonlat(lon, lat);
      Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

      int x = (int)DRG_pix[0];
      int y = (int)DRG_pix[1];

      PixelMask<PixelGray<uint8> > DRGVal = interpDRG(x, y);
      //insert data
      allImgPts[i] = (float) DRGVal;
    }else {
      //write -1 to designate and invalid point
      allImgPts[i] = -1;
    }
  }
  return allImgPts;
}

int main( int argc, char *argv[] ) {


  GlobalParams globalParams;
  //globalParams.reflectanceType = NO_REFL;
  globalParams.reflectanceType = LUNAR_LAMBERT;
  //globalParams.reflectanceType = LAMBERT;
  globalParams.slopeType = 1;
  globalParams.shadowThresh = 40;
  globalParams.albedoInitType = 1;
  globalParams.exposureInitType = 1;

  ModelParams modelParams;
  modelParams.exposureTime = 1.0;
  modelParams.rescalingParams[0] = 1;
  modelParams.rescalingParams[1] = 0;

  modelParams.sunPosition[0] = 1000*72987509.682619;//*sunPositions[i][0];
  modelParams.sunPosition[1] = 1000*133319340.07726;//*sunPositions[i][1];
  modelParams.sunPosition[2] = 1000*555820.93952321;//*sunPositions[i][2];

  modelParams.spacecraftPosition[0] = 1000*1668.1656423675;
  modelParams.spacecraftPosition[1] = 1000*127.53606522819;
  modelParams.spacecraftPosition[2] = 1000*774.17340580747;


  string inputCSVFilename; 
  string inputDEMFilename;
  string DRGFilename;  
  string DEMFilename; 

  bool aligned_set_up_tracks = false;
  bool update_model_params = true;


  aligned_set_up_tracks = false;
  inputCSVFilename = string("../data/Apollo15-LOLA/RDR_2E4E_25N27NPointPerRow_csv_table.csv"); 
  inputDEMFilename = string("../data/Apollo15-DEM/1134_1135-DEM.tif");
  
  //DRGFilename = string("../data/Apollo15-DRG/1134_1135-DRG.tif");  
  DRGFilename = string("../data/Apollo15-DRG/AS15-M-1134_map.tif");  
  DEMFilename = string("../results/dem.tiff"); 


  // Need this! Commented out for debugging! 
  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);
  /*printf("numTracks = %d\n", trackPts.size());
    for(int i = 0; i < trackPts.size(); i++){
    printf("numShots[%d] = %d\n", i, trackPts[i].size());
    }
   */

  int numVerPts = 6000;
  int numHorPts = 6000;
  vector<int> trackIndices;
  trackIndices.resize(trackPts.size());
  for (int i = 0; i < trackPts.size(); i++){
    trackIndices[i] = i;
  }
  // see need this above
  //write all tracks to image
  // Need this! Commented out for debugging! 
  //MakeGrid(trackPts, numVerPts, numHorPts, DEMFilename, trackIndices);

  //debugg GetTrackPtsFromImage by calling a little test function that looks very similar and records/prints out some important statistics.
  //WhenBuildImgPts(trackPts[1]);
  float normalizer_top = 0.0;
  float num_valid = 0.0;

  if( aligned_set_up_tracks ){
    for (int k = 1; k < trackPts.size(); k++){

    //string outFilename = "demo.tiff";
    //if (k == 1){
    //ShowTrackPtsOnImage(trackPts[k], DRGFilename, outFilename);
    //}
      std::string filename;
      char* filename_char = new char[500];
      sprintf (filename_char, "../results/dem_orbit_%d.txt", k);
      filename = std::string(filename_char);
      vector<float> pts = GetTrackPtsByID(trackPts[k], 3);
      SaveVectorToFile(pts, filename);
      
      /*
      vector<float> reflectance = ComputeTrackReflectance(trackPts[k], modelParams, globalParams);
      std::string reflectanceFilename;
      char* reflectanceFilename_char = new char[500];
      sprintf (reflectanceFilename_char, "../results/refl_orbit_%d.txt", k);
      reflectanceFilename = std::string(reflectanceFilename_char);
      SaveVectorToFile(reflectance, reflectanceFilename);
      printf("k = %d, reflectance.size() = %d\n", k, reflectance.size());
      //normalizer = normalizer_for_this_track(img_vals,refl_array);
      */
      /*
      vector<float> imgPts;
      imgPts = GetTrackPtsFromImage(trackPts[k], DRGFilename, 3);    
      std::string imgPtsFilename;
      char* imgPtsFilename_char = new char[500];
      sprintf (imgPtsFilename_char, "../results/img_orbit_%d.txt", k);
      imgPtsFilename = std::string(imgPtsFilename_char);
      SaveVectorToFile(imgPts, imgPtsFilename);
      */
      // We have two choices output: output image at every pt or every valid refl
      // allImagPts.size() = trackPts.size(), at points where the track is invalid we havea -1
      vector<float> allImgPts;
      allImgPts = AllTrackPtsFromImage(trackPts[k],DRGFilename);
      std::string allImgPtsFilename;
      char* allImgPtsFilename_char = new char[500]; 
      sprintf(allImgPtsFilename_char, "../results/all_img_track_orbit_%d.txt",k);
      allImgPtsFilename = std::string(allImgPtsFilename_char);
      SaveVectorToFile(allImgPts, allImgPtsFilename);

      cout << "Starting normalization calc." << endl;
      /*
      //account for normalization here
      for(int m = 0; m < reflectance.size(); m ++){
        if(reflectance[m]!= -1){
          normalizer_top += allImgPts[m];
          //normalizer_top += allImgPts[m]/reflectance[m];
          num_valid += 1.0;
        }
      }
      */
      float scaleFactor;
      /*
      //determine the scalling factor and save the synthetic image points- START
      
      scaleFactor = ComputeScaleFactor(allImgPts, reflectance);
      vector<float> synthImg = ComputeSyntImgPts(scaleFactor, reflectance);
      float matchingErr = ComputeMatchingError(synthImg, allImgPts);
      printf("matchingError  %f\n", matchingErr);
      
      std::string synthImgPtsFilename;
      char* synthImgPtsFilename_char = new char[500];
      sprintf (synthImgPtsFilename_char, "../results/synthImg_track_%d.txt", k);
      synthImgPtsFilename = std::string(synthImgPtsFilename_char);
      SaveVectorToFile(synthImg, synthImgPtsFilename);
      //determine the scalling factor - END
      */
      printf("\n\nInfo: k = %d, normalizer_top = %f, num_valid = %f, scaleFactor = %f\n",
          k, normalizer_top, num_valid, scaleFactor);
      //reflectance.clear();
      //
      vector<float> demPts;
      demPts = GetTrackPtsFromDEM(trackPts[k], inputDEMFilename, 3);
      std::string demPtsFilename;
      char* demPtsFilename_char = new char[500];
      sprintf (demPtsFilename_char, "../results/sdem_orbit_%d.txt", k);
      demPtsFilename = std::string(demPtsFilename_char);
      SaveVectorToFile(demPts, demPtsFilename);

      //write individual tracks to image
      string outDEMFilename;
      char* outDEMFilename_char = new char[500];
      sprintf (outDEMFilename_char, "../results/dem_%d.tiff", k);
      outDEMFilename = std::string(outDEMFilename_char);
      printf("demfilename = %s\n", outDEMFilename.c_str());
      vector<int> trackIndices;
      trackIndices.resize(1);
      trackIndices[0] = k;
      MakeGrid(trackPts, numVerPts, numHorPts, outDEMFilename, trackIndices);

      delete[] outDEMFilename_char;
      delete[] demPtsFilename_char;
      delete[] allImgPtsFilename_char;
      //delete[] imgPtsFilename_char;
      //delete[] reflectanceFilename_char;
      delete[] filename_char;

    }
  }

  if( update_model_params){
 
    int maxNumIter = 10;
    int maxNumStarts = 160;
    vector<Vector<float, 6> >init_d_array;
    init_d_array.resize(maxNumStarts);
    for (int i = 0; i < maxNumStarts; i++){
      init_d_array[i][0] = 1.0;
      init_d_array[i][1] = 0.0;
      init_d_array[i][2] = (i-maxNumStarts/2)*5;
      init_d_array[i][3] = 0.0;
      init_d_array[i][4] = 1.0;
      init_d_array[i][5] = 0.0;//(i-maxNumStarts/2)*25;
    }    

    vector<Vector<float, 6> > final_d_array;
    vector<float> error_array;
    error_array.resize(maxNumStarts);
    final_d_array.resize(maxNumStarts);
  
    
    DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
    GeoReference DRGGeo;
    read_georeference(DRGGeo, DRGFilename);
  

    ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
                                                                         ConstantEdgeExtension()),
  							    BilinearInterpolation());

    cout <<"Done interpolating the subsampled image"<<endl;
  

    //get the true image points
    cout << "GetAllPtsFromImage..." << endl; 
    GetAllPtsFromImage(trackPts, interpDRG, DRGGeo);
   
    cout << "ComputeTrackReflectance..." << endl;
    ComputeAllReflectance(trackPts, modelParams, globalParams);

    cout << "Calling UpdateMatchingParams ..."<< endl;
    //return matching error and transform
    UpdateMatchingParams(trackPts, DRGFilename, 
                         modelParams, globalParams,maxNumIter,  
                         init_d_array, final_d_array, error_array);

    int bestResult = 0;
    float smallestError = error_array[0];
    for (int index = 0; index < init_d_array.size(); index++){
        printf("OUT %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
	      index, error_array[index], 
              final_d_array[index](0), final_d_array[index](1), 
              final_d_array[index](2), final_d_array[index](3),
	      final_d_array[index](4), final_d_array[index](5));
        if  (error_array[index] < smallestError){
	     smallestError = error_array[index]; 
             bestResult = index;
        }      
    }    

    //write results to image outside matching
    string DRGFilenameNoPath = sufix_from_filename(DRGFilename);
    string outFilename = "../results" + prefix_less3_from_filename(DRGFilenameNoPath) + "_results.tif";  
    printf("bestResult = %d\n", bestResult);
    ShowFinalTrackPtsOnImage(trackPts, final_d_array[bestResult], 
                             trackIndices, DRGFilename, outFilename);
    cout << "UpdateMatchingParams finshed." << endl;
  }
  /*
  //save the, # valid, normalizer, & the division
  std::string matchingParamsFilename;
  char* matchingParamsFilename_char = new char[500];
  sprintf( matchingParamsFilename_char, "../results/matching_params");
  matchingParamsFilename = std::string(matchingParamsFilename);

  // the following is very silly
  vector<float> read_out;
  for(int i=0; i< d.size(); i++ )
  {
  read_out.push_back(d[i]);
  }
  SaveVectorToFile(read_out,matchingParamsFilename);
   */

  /*
     vector<float> norm_consts;
     norm_consts.push_back( normalizer_top);
     norm_consts.push_back(  num_valid);
     norm_consts.push_back( normalizer_top/num_valid);
     std::string normConstsFilename;
     char* normConstsFilename_char = new char[500];
  //sprintf (normConstsFilename_char, "../results/save_normalizer_constants.txt");
  normConstsFilename = std::string(normConstsFilename_char);
  SaveVectorToFile(norm_consts, normConstsFilename);

  delete normConstsFilename_char;
   */

  //validation...?
  // 'basin of attraction'...
  // how sensitive is the solution to the initialization
  // i.e. if we shift the initial position (- 10, + 5) pixels does this affect the final solution.  Can we move +/- 2000 pixels and still get the same solution?  I.e. is the problem convex.

  // working to accomplish this by changing the update params interface



  return 0;
}


















