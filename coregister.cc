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

#include "coregister.h"

int GetTimeDiff(pointCloud prevPt, pointCloud currPt, float timeThresh)
{


  double prevTime = prevPt.hour*3600.0 + prevPt.min*60.0 + prevPt.sec;
  double currTime = currPt.hour*3600.0 + currPt.min*60.0 + currPt.sec;
 
  if (prevPt.year != currPt.year){//new orbit
      return (1);
  }  

  if (prevPt.month != currPt.month){ //new orbit
      return (1);
  } 
  
  if (prevPt.day != currPt.day){//new orbit
      return (1);
  } 
 
  if (currTime - prevTime > timeThresh){ //new orbit
     return (1);
  }
  else{
     return (0);
  }
}
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


void SaveVectorToFile(vector<float> v, string filename)
{
 

  FILE *fp;

  fp = fopen(filename.c_str(), "w");

  for (int i = 0; i < v.size()-1; i++){
    fprintf(fp, "%f\n", v[i]);
  }
  fprintf(fp, "%f", v[v.size()-1]);
  fclose(fp);
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
vector<float> GetOrbitPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename, int ID)
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

vector<float> GetOrbitPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID)
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


void UpdateMatchingParams(string DRGFilename)
{

   DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
   ImageView<float> x_deriv = derivative_filter(DRG, 1, 0);
   ImageView<float> y_deriv = derivative_filter(DRG, 0, 1);
   
   Vector<float,6> d;//defines the affine transform
   d(0) = 1.0;
   d(4) = 1.0;
   Matrix<float,6,6> rhs;
   Vector<float,6> lhs;

   int ii, jj, x_base, y_base;
   float xx = x_base + d[0] * ii + d[1] * jj + d[2];
   float yy = y_base + d[3] * ii + d[4] * jj + d[5];
   
   // Fill in symmetric entries
   rhs(1,0) = rhs(0,1);
   rhs(2,0) = rhs(0,2);
   rhs(2,1) = rhs(1,2);
   rhs(1,3) = rhs(0,4);
   rhs(2,3) = rhs(0,5);
   rhs(2,4) = rhs(1,5);
   rhs(3,0) = rhs(0,3);
   rhs(3,1) = rhs(1,3);
   rhs(3,2) = rhs(2,3);
   rhs(4,0) = rhs(0,4);
   rhs(4,1) = rhs(1,4);
   rhs(4,2) = rhs(2,4);
   rhs(4,3) = rhs(3,4);
   rhs(5,0) = rhs(0,5);
   rhs(5,1) = rhs(1,5);
   rhs(5,2) = rhs(2,5);
   rhs(5,3) = rhs(3,5);
   rhs(5,4) = rhs(4,5);
   try {
       solve_symmetric_nocopy(rhs,lhs);
   } catch (ArgumentErr &/*e*/) {
                //             std::cout << "Error @ " << x << " " << y << "\n";
                //             std::cout << "Exception caught: " << e.what() << "\n";
                //             std::cout << "PRERHS: " << pre_rhs << "\n";
                //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
                //             std::cout << "RHS: " << rhs << "\n";
                //             std::cout << "LHS: " << lhs << "\n\n";
                //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
                //             exit(0);
   }
   d += lhs;
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

vector<vector<LOLAShot> > CSVFileRead(string CSVFilename)
{
  string line;
  ifstream myfile (CSVFilename.c_str());
  int lineIndex = 0;
  //int shotIndex = 0;
 
  int trackIndex;
  vector<vector<LOLAShot> >trackPts;
  LOLAShot shot;
 
  pointCloud currPt;
  pointCloud prevPt;

  if (myfile.is_open())
  {
    while (! myfile.eof() )
    {
      getline (myfile,line);
     
      if (lineIndex > 0){//skip header
     
        //float lon, lat, rad;
        char *temp = new char[160]; 
        char *lon = new char[160]; 
        char *lat = new char[160]; 
        char *rad = new char[160];
        char *detID = new char[5];
        char *tmpf = new char[200];
	
        //printf("line = %s\n", line.c_str()); 
	sscanf(line.c_str(), "%s %s %s %s %s %s %s %s %s %s  %s", 
	       temp, lon, lat, rad, tmpf, tmpf, tmpf, tmpf, tmpf, tmpf, detID);

        string time = temp;
        char year[5]; 
        char month[3]; 
        char day[3];
        char hour[3];
        char min[3];
        char sec[12];
        char s[2];
	
	string detIDs = detID;

        if (time.length() > 1){
	  size_t length;
	  length = time.copy(year, 4, 0);
          //printf("length = %d\n", length);
	  year[length] = '\0';
          currPt.year = atof(year); 
       
	  length = time.copy(month, 2, 5);
          //printf("length = %d\n", length);
	  month[length] = '\0';
          currPt.month = atof(month); 

	  length = time.copy(day, 2, 8);
	  day[length] = '\0';  
          currPt.day = atof(day);

	  length = time.copy(hour, 2, 11);
	  hour[length] = '\0';
          currPt.hour = atof(hour);

	  length = time.copy(min, 2, 14);
	  min[length] = '\0';
	  length = time.copy(sec, 11, 17);
	  sec[length] = '\0';
          currPt.sec = atof(sec);

          length = detIDs.copy(s, 1, 2);
          s[length] = '\0';
          currPt.s = atof(s);
	  //printf("%s %s %s %s %s %s detID = %s\n", year, month, day, hour, min, sec, s);
	}
	
        Vector3 coords;         
        currPt.coords(0) = atof(lon);
        currPt.coords(1) = atof(lat);
	currPt.coords(2) = atof(rad);
        
        
        if ((currPt.coords(0)!=0.0) && (currPt.coords(1)!=0.0) ){ //valid lidar point
     
	  if (lineIndex == 1){ //initialize first track
	      trackIndex = 0;
	      trackPts.resize(trackIndex+1);
              printf("lineIndex = %d\n", lineIndex);
          }
          else{
            
	     if (GetTimeDiff(prevPt, currPt, 3000)){ //new track
	         trackPts[trackIndex].push_back(shot);//add last shot to the previous track
                 shot.LOLAPt.clear();
                 trackIndex++;
 	         trackPts.resize(trackIndex+1); //start new track
                 shot.LOLAPt.push_back(currPt);
	     }
	     else{ //same track
                 if (GetTimeDiff(prevPt, currPt, 0)){//new shot
	             trackPts[trackIndex].push_back(shot);
                     shot.LOLAPt.clear();
                     shot.LOLAPt.push_back(currPt);
                 }
	         else{ //same shot
                    shot.LOLAPt.push_back(currPt);
                 }
             }
	   }

      
           //copy current pc into prevPt
           prevPt.coords(0) = currPt.coords(0);
           prevPt.coords(1) = currPt.coords(1);
           prevPt.coords(2) = currPt.coords(2);
           prevPt.year = currPt.year;
           prevPt.month = currPt.month;
           prevPt.day = currPt.day;
           prevPt.hour = currPt.hour;
           prevPt.min = currPt.min;
           prevPt.sec = currPt.sec;   
           prevPt.s = currPt.s;
	}

        delete temp;
        delete lon;
	delete lat;
	delete rad;
	delete detID;
      } 
      lineIndex++; 
    }
    myfile.close();
  }

  else cout << "Unable to open file";

  
  return trackPts; 
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


  //#if 0
  
  string inputCSVFilename = string("../data/Apollo15-LOLA/RDR_2E4E_25N27NPointPerRow_csv_table.csv"); 
  string inputDEMFilename = string("../data/Apollo15-DEM/1134_1135-DEM.tif");
  string DRGFilename = string("../data/Apollo15-DRG/1134_1135-DRG.tif");  
  string DEMFilename = string("../results/dem.tiff"); 

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

  for (int k = 0; k < trackPts.size(); k++){
    

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
    reflectance.clear();

    vector<float> imgPts;
    imgPts = GetOrbitPtsFromImage(trackPts[k], DRGFilename, 3);    
    std::string imgPtsFilename;
    char* imgPtsFilename_char = new char[500];
    sprintf (imgPtsFilename_char, "../results/img_orbit_%d.txt", k);
    imgPtsFilename = std::string(imgPtsFilename_char);
    //SaveTrackImgPts(imgPts, imgPtsFilename);
     SaveVectorToFile(imgPts, imgPtsFilename);

    vector<float> demPts;
    demPts = GetOrbitPtsFromDEM(trackPts[k], inputDEMFilename, 3);
    std::string demPtsFilename;
    char* demPtsFilename_char = new char[500];
    sprintf (demPtsFilename_char, "../results/sdem_orbit_%d.txt", k);
    demPtsFilename = std::string(demPtsFilename_char);
    //SaveTrackImgPts(demPts, demPtsFilename);
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
    
  }


  //#endif

  #if 0
  string inputDEMFile = string("/Users/anefian/projects/stereopipeline/data/lola/ldem_45N_100m.tif");
  string synthImgFile = string("../results/synth_45N_100m.tif");
  DiskImageView<PixelGray<uint8> >   inputDEM(inputDEMFile); 
  GeoReference inputDEMGeo;
  read_georeference(inputDEMGeo, inputDEMFile);

  //ImageView<PixelMask<PixelGray<float> > > synthImg (inputDEM.cols(), inputDEM.rows());  

  //TODO: determine the bounding box (XYZ coordinates) for the lola data

  //TODO: read the real image

  //TODO: create an artificial image
  for (int k = 1 ; k < 1000/*inputDEM.rows()*/; k++) {
    for (int l = 1; l < 1000/*inputDEM.cols()*/; l++) {     

	 printf("l = %d, k = %d\n", l, k); 
	 Vector2 lonlat2; 
	 //get the XYZ coordinates of the current point        
         Vector2 inputDEMVal(l,k);
         lonlat2 = inputDEMGeo.pixel_to_lonlat(inputDEMVal);
         Vector3 lonlat3(lonlat2(0),lonlat2(1),(inputDEM)(l,k));
	 Vector3 xyz = inputDEMGeo.datum().geodetic_to_cartesian(lonlat3);
         printf("xyz[0]=%f, xyz[1]=%f, xyz[2]=%f\n", xyz(0), xyz(1), xyz(2));

         //get the XYZ coordinates of the left point  
	 Vector2 inputDEMLeftVal;
	 inputDEMLeftVal(0) = l-1;
	 inputDEMLeftVal(1) = k;
	 lonlat2 = inputDEMGeo.pixel_to_lonlat(inputDEMLeftVal);
	 Vector3 lonlat3Left(lonlat2(0),lonlat2(1),(inputDEM)(l-1,k));
	 Vector3 xyzLeft = inputDEMGeo.datum().geodetic_to_cartesian(lonlat3Left);
	 
	 //get the XYZ coordinates of the top point  
	 Vector2 inputDEMTopVal;
	 inputDEMTopVal(0) = l;
	 inputDEMTopVal(1) = k-1;
         lonlat2 = inputDEMGeo.pixel_to_lonlat(inputDEMTopVal);
         Vector3 lonlat3Top(lonlat2(0),lonlat2(1),(inputDEM)(l,k-1));
	 Vector3 xyzTop = inputDEMGeo.datum().geodetic_to_cartesian(lonlat3Top);

         Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);
         float reflectance = ComputeReflectance(normal, xyz, t_modelParams, globalParams);
         printf("n[0]=%f, n[1]=%f, n[2]=%f, reflectance = %f\n", normal(0), normal(1), normal(2), reflectance);
      
         //synthImg(l,k) = reflectance;
      }
  }
  #endif
  /*
  //write the synthetic image
  write_georeferenced_image(synthImgFile, 
                            channel_cast<uint8>(clamp(synthImg,0.0,255.0)),
                            inputDEMGeo, TerminalProgressCallback());
  */
  //TODO: determine the feature points in the synthetic image
  
  //TODO: determine the feature points in the real image
}


















