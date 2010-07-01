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
//#include <cv.h>
//#include <highgui.h>
#include "io.h"
#include "match.h"
#include "coregister.h"


vector<float> GetTrackPtsByID(vector<LOLAShot> trackPts, int ID)
{
  vector<float> pts;
  for (int i = 0; i < trackPts.size(); i++){
    for (int k = 0; k < trackPts[i].LOLAPt.size(); k++){

      float rad = trackPts[i].LOLAPt[k].coords[2];
      int id = trackPts[i].LOLAPt[k].s;

      if (id == ID){
        pts.push_back(rad);
      }

    }
  }

  return pts;
}

Vector4 FindMinMaxLat(vector<vector<LOLAShot> >trackPts)
{
  float minLat = 180;
  float maxLat = -180;
  float minLon = 180;
  float maxLon = -180;

  for (int i = 0; i < trackPts.size(); i++){
    for (int j = 0; j < trackPts[i].size(); j++){
      for(int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){
        float lon = trackPts[i][j].LOLAPt[k].coords[0];
        float lat = trackPts[i][j].LOLAPt[k].coords[1];

        if (lat < minLat){
          minLat = lat;
        }

        if (lon < minLon){
          minLon = lon;
        }

        if (lat > maxLat){
          maxLat = lat;
        }

        if (lon > maxLon){
          maxLon = lon;
        }
      }
    }
  }

  Vector4 coords;
  coords(0) = minLat; 
  coords(1) = maxLat; 
  coords(2) = minLon; 
  coords(3) = maxLon;

  return coords;

}

pointCloud GetPointFromIndex(vector<pointCloud> LOLAPts, int index)
{
  pointCloud pt;
  pt.s = -1;//invalid pointCloud
  for(int i = 0;i < LOLAPts.size(); i++){
    if (LOLAPts[i].s == index){
      return LOLAPts[i];
    }
  }

  return pt;
}
Vector3 ComputeNormal(vector<pointCloud> LOLAPts)
{
  Vector3 normal;
  Matrix<float,5,3> A;
  for(int i = 0;i < LOLAPts.size(); i++){
    //transform into x,y,z coordinates
    GeoReference DEMGeo;
    Vector3 xyz = DEMGeo.datum().geodetic_to_cartesian(LOLAPts[i].coords);

    A(i,0) = xyz[0];
    A(i,1) = xyz[1];
    A(i,2) = xyz[2];
  }
  return normal;
}

float ComputeReflectance(vector<pointCloud> LOLAPts, ModelParams modelParams, GlobalParams globalParams)
{
  pointCloud centerPt = GetPointFromIndex(LOLAPts, 3);
  pointCloud topPt = GetPointFromIndex(LOLAPts, 2);
  pointCloud leftPt = GetPointFromIndex(LOLAPts, 1);

  if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){
    Datum moon;
    moon.set_well_known_datum("D_MOON");

    centerPt.coords[2] = (centerPt.coords[2]-1737.4)*1000;
    topPt.coords[2] = (topPt.coords[2]-1737.4)*1000;
    leftPt.coords[2] = (leftPt.coords[2]-1737.4)*1000;
    printf("c = %f, t = %f, l = %f\n", centerPt.coords[2],topPt.coords[2],leftPt.coords[2]);

    Vector3 xyz = moon.geodetic_to_cartesian(centerPt.coords);
    Vector3 xyzTop = moon.geodetic_to_cartesian(topPt.coords);
    Vector3 xyzLeft = moon.geodetic_to_cartesian(leftPt.coords);
    printf("xyz = %f %f %f\n", xyz(0), xyz(1), xyz(2));
    Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);
    //printf("normal = %f %f %f\n", normal(0), normal(1), normal(2));
    float reflectance = ComputeReflectance(normal, xyz, modelParams, globalParams);

    return reflectance;
  }
  else{
    return -1;
  }
}
vector<float> ComputeTrackReflectance(vector<LOLAShot> trackPts, ModelParams modelParams, GlobalParams globalParams)
{
  vector<float> reflectance;
  reflectance.resize(trackPts.size());
  for (int m = 0; m < trackPts.size();m++){
    reflectance[m] = ComputeReflectance(trackPts[m].LOLAPt, modelParams, globalParams);
    printf("ref = %f\n", reflectance[m]);
  }
  return reflectance;
}

void WhenBuildImgPts(vector<LOLAShot> trackPts){
 //this function is to understand which pts GetTrackPtsFromImage will extract
 
 /* First set of tests:
    1. what is the length of trackPts[0] // [0] is the index being fed in
    2. what is the # accepted based on the id == 3 condition
    3. what is the # accepted on the valid refl. condition?
 */
 
  //1. # of track points
  printf("\n\n# trackPts = %d\n",trackPts.size());
  
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


vector<float> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename, int ID)
{
  DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  vector<float> imgPts;
  //imgPts.resize(orbitPts.size());

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  for (int i = 0; i < trackPts.size(); i++){
    for (int k = 0; k < trackPts[i].LOLAPt.size(); k++){

      float lon = trackPts[i].LOLAPt[k].coords[0];
      float lat = trackPts[i].LOLAPt[k].coords[1];
      float rad = trackPts[i].LOLAPt[k].coords[2];
      int id = trackPts[i].LOLAPt[k].s;

      Vector2 DEM_lonlat(lon, lat);
      Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

      int x = (int)DRG_pix[0];
      int y = (int)DRG_pix[1];

      PixelMask<PixelGray<uint8> > DRGVal = interpDRG(x, y);
      if (id == ID){
        imgPts.push_back((float)DRGVal);
      }

    }
  }

  return imgPts;
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


vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID)
{
  DiskImageView<PixelGray<float> >   DEM(DEMFilename);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, DEMFilename);

  vector<float> demPts;

  ImageViewRef<PixelGray<float>  >  interpDEM = interpolate(edge_extend(DEM.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  for (int i = 0; i < trackPts.size(); i++){
    for(int j = 0; j < trackPts[i].LOLAPt.size(); j++){
      float lon = trackPts[i].LOLAPt[j].coords[0];
      float lat = trackPts[i].LOLAPt[j].coords[1];
      float rad = trackPts[i].LOLAPt[j].coords[2];
      int id = trackPts[i].LOLAPt[j].s;

      Vector2 DEM_lonlat(lon, lat);
      Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(DEM_lonlat);

      int x = (int)DEM_pix[0];
      int y = (int)DEM_pix[1];

      PixelGray<float>  DEMVal = interpDEM(x, y);

      if (id == ID){
        demPts.push_back((float)DEMVal);
      }             
    }
  }

  return demPts;
}

void MakeGrid(vector<vector<LOLAShot> >trackPts, int numVerPts, int numHorPts, string DEMFilename, vector<int> trackIndices)
{
  int l, m, n;  
  ImageView<PixelGray<float> > DEMImage(numHorPts, numVerPts);
  GeoReference DEMgeo;

  Vector4 coords = FindMinMaxLat(trackPts);

  printf("minLat=%f, maxLat=%f, minLon=%f maxLon=%f\n", coords(0), coords(1), coords(2), coords(3));

  float minLat = coords(0); 
  float maxLat = coords(1); 
  float minLon = coords(2); 
  float maxLon = coords(3);

  float lonDelta = (maxLon-minLon)/numHorPts;
  float latDelta = (maxLat-minLat)/numVerPts;

  //printf("lonDelta = %f, latDelta = %f\n", lonDelta, latDelta);
  //init the DEM
  for (l = 0; l < numHorPts; l++){
    for (m = 0; m < numVerPts; m++){
      DEMImage(l, m) = 0.0;
    }
  }
  //fill the DEM

  printf("numTracks = %d\n", trackIndices.size());
  for (int k = 0; k < trackIndices.size();k++){
    int trackIndex = trackIndices[k];
    printf("trackIndex = %d\n", trackIndex);
    for (n = 0; n < trackPts[trackIndex].size(); n++){ 
      for (int s = 0; s < trackPts[trackIndex][n].LOLAPt.size(); s++){

        float lon_index = (trackPts[trackIndex][n].LOLAPt[s].coords[0] - minLon)/lonDelta;
        float lat_index = (trackPts[trackIndex][n].LOLAPt[s].coords[1] - minLat)/latDelta;
        l = (int)floor(lon_index);
        m = (int)floor(lat_index);

        if ((m < numVerPts) && (l<numHorPts)){ 
          DEMImage(l, m) = trackPts[trackIndex][n].LOLAPt[s].coords[2]; 
        }
        else{
          printf("Error\n");
          printf("l = %d, m = %d, numHorPts = %d, numVerPts = %d\n", l, m, numHorPts, numVerPts);
        }
      }
    }
  }
  write_georeferenced_image(DEMFilename, 
      DEMImage,
      DEMgeo, TerminalProgressCallback("{Core}","Processing:"));

}

float normalizer_for_this_track(vector<float> img_vals,vector<float> refl_array){
  float n_here = 0;

  for(int i = 0; i < refl_array.size(); i++){
    if(refl_array[i]!= -1){
      printf("If statement in: norm... works");
      n_here += img_vals[i];
     
    }
  }

  return n_here;
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

  /*
  int comp_number = 1; // 0 = Ara's paths, 1 = Dave's paths
  if (comp_number == 0)
  {
  string inputCSVFilename = string("../data/Apollo15-LOLA/RDR_2E4E_25N27NPointPerRow_csv_table.csv"); 
  string inputDEMFilename = string("../data/Apollo15-DEM/1134_1135-DEM.tif");
  string DRGFilename = string("../data/Apollo15-DRG/1134_1135-DRG.tif");  
  string DEMFilename = string("../results/dem.tiff"); 
  }
  if(comp_number == 1)
  {
  */
  string inputCSVFilename = string("../../data/Apollo15-LOLA/RDR_2E4E_25N27NPointPerRow_csv_table.csv"); 
  string inputDEMFilename = string("../../data/Apollo15-DEM/1134_1135-DEM.tif");
  string DRGFilename = string("../../data/Apollo15-DRG/1134_1135-DRG.tif");  
  string DEMFilename = string("../../results/dem.tiff"); 
  //}
  
  vector<vector<LOLAShot> > trackPts =  CSVFileRead(inputCSVFilename);
  printf("numTracks = %d\n", trackPts.size());
  for(int i = 0; i < trackPts.size(); i++){
    printf("numShots[%d] = %d\n", i, trackPts[i].size());
  }

  int numVerPts = 6000;
  int numHorPts = 6000;
  vector<int> trackIndices;
  trackIndices.resize(trackPts.size());
  for (int i = 0; i < trackPts.size(); i++){
    trackIndices[i] = i;
  }

  //write all tracks to image
  MakeGrid(trackPts, numVerPts, numHorPts, DEMFilename, trackIndices);

  //debugg GetTrackPtsFromImage by calling a little test function that looks very similar and records/prints out some important statistics.
  WhenBuildImgPts(trackPts[1]);
  float normalizer_top = 0.0;
  float num_valid = 0.0;
  for (int k = 1; k < trackPts.size(); k++){


    std::string filename;
    char* filename_char = new char[500];
    sprintf (filename_char, "../results/dem_orbit_%d.txt", k);
    filename = std::string(filename_char);
    vector<float> pts = GetTrackPtsByID(trackPts[k], 3);
    SaveVectorToFile(pts, filename);

    vector<float> reflectance = ComputeTrackReflectance(trackPts[k], modelParams, globalParams);
    std::string reflectanceFilename;
    char* reflectanceFilename_char = new char[500];
    sprintf (reflectanceFilename_char, "../results/refl_orbit_%d.txt", k);
    reflectanceFilename = std::string(reflectanceFilename_char);
    SaveVectorToFile(reflectance, reflectanceFilename);
    printf("k = %d, reflectance.size() = %d\n", k, reflectance.size());
    //normalizer = normalizer_for_this_track(img_vals,refl_array);
      
    vector<float> imgPts;
    imgPts = GetTrackPtsFromImage(trackPts[k], DRGFilename, 3);    
    std::string imgPtsFilename;
    char* imgPtsFilename_char = new char[500];
    sprintf (imgPtsFilename_char, "../results/img_orbit_%d.txt", k);
    imgPtsFilename = std::string(imgPtsFilename_char);
    SaveVectorToFile(imgPts, imgPtsFilename);

    // We have two choices output: output image at every pt or every valid refl
    // allImagPts.size() = trackPts.size(), at points where the track is invalid we havea -1
    vector<float> allImgPts;
    allImgPts = AllTrackPtsFromImage(trackPts[k],DRGFilename);
    std::string allImgPtsFilename;
    char* allImgPtsFilename_char = new char[500]; 
    sprintf(allImgPtsFilename_char, "../results/all_img_track_orbit_%d.txt",k);
    allImgPtsFilename = std::string(allImgPtsFilename_char);
    SaveVectorToFile(allImgPts, allImgPtsFilename);

    //account for normalization here
    for(int m = 0; m < reflectance.size(); m ++){
      if(reflectance[m]!= -1){
        normalizer_top += allImgPts[m];
        num_valid += 1.0;
      }
    }
    printf("\n\nInfo: k = %d, normalizer_top = %f, num_valid = %f\n",k,normalizer_top,num_valid);
    //reflectance.clear();
    
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
    delete[] imgPtsFilename_char;
    delete[] reflectanceFilename_char;
    delete[] filename_char;
    
  }

 // UpdateMatchingParams(trackPts, DRGFilename, modelParams, globalParams);

  //save the, # valid, normalizer, & the division
  
  vector<float> norm_consts;
  norm_consts.push_back( normalizer_top);
  norm_consts.push_back(  num_valid);
  norm_consts.push_back( normalizer_top/num_valid);
  std::string normConstsFilename;
  char* normConstsFilename_char = new char[500];
  sprintf (normConstsFilename_char, "../results/save_normalizer_constants.txt");
  normConstsFilename = std::string(normConstsFilename_char);
  SaveVectorToFile(norm_consts, normConstsFilename);
 
  delete[] normConstsFilename_char;
}


















