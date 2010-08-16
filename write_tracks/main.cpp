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

#include "main.h"

int GetTimeDiff(pointCloud prevPt, pointCloud currPt, float timeThresh)
{
  //dtj June 18th 2010 - reading code.  So we're defining an 'orbit' as a single day?

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

void Min_Max_Pixel(int point_size, int &pix_R,int &pix_C,int &min_R,int &max_R,int &min_C,int &max_C){
  if( (pix_R - point_size) < min_R){
    min_R = pix_R - point_size;
  }
  if( (pix_R + point_size) > max_R){
    max_R = pix_R + point_size;
  }
  if( (pix_C - point_size) < min_C){
    min_C = pix_C - point_size;
  }
  if( (pix_C + point_size) > max_C){
    max_C = pix_C + point_size;
  }
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
/*
Vector4 Find_crop(vector<vector<LOLAShot > > trackPts, string DRGName ){
    
    // internal variables
    DRG  

    // load up drg reference - ideally pass by reference

    // calc. min/max lat & lon

    // translate coordinates to pixels

    // save the pixels
}
*/
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
      if(lineIndex==46821){
        printf("Danger: lineIndex = %d",lineIndex);
        cout << "Get warning off before?";
      }
      getline (myfile,line);

      if(!myfile.eof()){




        if(lineIndex==46821){
          printf("We didn't get this far - did we?");
        }     

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
            //printf("Copy Current point, lineIndex = %d\n",lineIndex);
          }

          delete temp;
          delete lon;
          delete lat;
          delete rad;
          delete detID;
        } 
        lineIndex++; 
        //printf("Copy Current Point, idx = %d\n",lineIndex);
      }
    }
    printf("Finished copying one files worth of points, lineIndex = %d\n",lineIndex);
    myfile.close();
  }

  else cout << "Unable to open file";

  printf("Finished copying points");  
  return trackPts; 
}

int main( int argc, char *argv[] ) {
/*
	//
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

  modelParams.sunPosition[0] = 1000*72987509.682619;//$*sunPositions[i][0];
  modelParams.sunPosition[1] = 1000*133319340.07726;//$*sunPositions[i][1];
  modelParams.sunPosition[2] = 1000*555820.93952321;//$*sunPositions[i][2];

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
  //MakeGrid(trackPts, numVerPts, numHorPts, DEMFilename, trackIndices);


*/

  //Overlay tracks
    //read DRG
    int data_source = 2;
    string drg_filename;
    string d_load_filename;
    if(data_source == 1){
      // string drg_filename = "1134_1135-DRG.tif";   
      //drg_filename = "../lima_svn_code/data/Apollo15-DRG/1134_1135-DRG.tif";
      //d_load_filename = "../lima_svn_code/results/1134_1135-d_final.txt";
    
      drg_filename = "../../data/Apollo15-DRG/1134_1135-DRG.tif";
      d_load_filename = "../../results/1134_1135-d_final.txt";
    
    }else if(data_source == 2){
      //drg_filename = "../lima_svn_code/data/Apollo15-DRG/AS15-M-1134_map.tif";
      //d_load_filename = "../lima_svn_code/results/AS15-M-1134_d_final.txt";
    
      drg_filename = "../../data/Apollo15-DRG/AS15-M-1134_map.tif";
      d_load_filename = "../../results/AS15-M-1134_d_final.txt";
    }
    

    //construct save names for these images
    std::string temp_name_sNow = sufix_from_filename(drg_filename);
    cout << "temp_name_sNow = " << temp_name_sNow << endl;;
    temp_name_sNow.erase(0,1);
    cout << "temp_name_sNow = " << temp_name_sNow << endl;; 
    string orginal_tracks =  prefix_less3_from_filename(temp_name_sNow) + "base_tracks.tif"; 
    string orginal_and_moved_tracks = prefix_less3_from_filename(temp_name_sNow) + "d_tracks.tif";  
 
    cout << "temp_name_sNow = " << temp_name_sNow << endl; 
    cout << "original_tracks = " << orginal_tracks << endl;
    cout << "original_and_moved_tracks = " << orginal_and_moved_tracks << endl;
                                           

    cout << "drg_filename = " << drg_filename << endl;
    cout << "d_load_filename = " << d_load_filename << endl;
             
    //load d file & print out
    FILE* d_final;
    float d0, d1, d2, d3, d4, d5;
    d_final = fopen(d_load_filename.c_str(),"r");    
    //d_final = fopen("AS-15_d_mod.txt","r");
    fscanf(d_final,"%*s %f %*s %f %*s %f %*s %f %*s %f %*s %f", &d0, &d1, &d2, &d3, &d4, &d5);
    printf("d[0] = %f d[1] = %f d[2] = %f\nd[3] = %f d[4] = %f d[5] = %f\n",d0,d1,d2,d3,d4,d5);
    //d_final.close();
    
    string csv_filename = "../data/Apollo15-LOLA/RDR_2E4E_25N27NPointPerRow_csv_table.csv";

    ImageView<PixelRGB<uint8> > drg_img;
    DiskImageResourceGDAL drg_rsrc(drg_filename);
    GeoReference drg_georef;

    read_georeference(drg_georef, drg_rsrc);
    read_image(drg_img, drg_filename);

    vector<vector<LOLAShot> > track_list = CSVFileRead(csv_filename);

    cout << endl;

    // orignal track only crop
    int min_R =1000000000 , max_R = 0 , min_C = 100000000, max_C = 0;
    int row = 0;
    int col = 0;

    int tracknum = 0;
    int numpix = 0;
    int point_size = 5;
    typedef vector<vector<LOLAShot> >::const_iterator track_itr_t;
    typedef vector<LOLAShot>::const_iterator shot_itr_t;
    typedef vector<pointCloud>::const_iterator point_itr_t;
    for (track_itr_t track = track_list.begin(); track != track_list.end(); track++) {
      for (shot_itr_t shot = track->begin(); shot != track->end(); shot++) {
        for (point_itr_t point = shot->LOLAPt.begin(); point != shot->LOLAPt.end(); point++) {
          Vector2 pix = drg_georef.lonlat_to_pixel(SubVector<const Vector3>(point->coords, 0, 2));
          if (BBox2(point_size, point_size, drg_img.cols() - point_size, drg_img.rows() - point_size).contains(pix)) {
            fill(crop(drg_img, int32(pix.x()) - point_size, int32(pix.y()) - point_size, point_size, point_size), PixelRGB<uint8>(0, 255, 255));
            //update current values for cropping
            row = int32( pix.y());
            col = int32( pix.x());
            Min_Max_Pixel( point_size, row, col, min_R, max_R, min_C, max_C);
            numpix++;
          }
        }
      }
      cout << "Status: " << ++tracknum << " / " << track_list.size() << "\r";
    }
    
    printf("Crop is rows ( %d, %d) cols( %d, %d)\n", min_R, max_R, min_C, max_C);
    cout << endl << "numpix = " << numpix << endl;
    /*
    for (out_itr track; 
      for (int j = 0;
    
    
    } */
    
    // write_image("1134_1135-base_tracks.tif", crop(drg_img, 1727, 1386, 3784, 4345), TerminalProgressCallback("vw", "BLAH: "));
    //write_image("/1134_1135-base_tracks.tif", crop(drg_img, 1727, 1386, 3784, 4345), TerminalProgressCallback("vw", "BLAH: "));
    //write_image(orginal_tracks.c_str(), crop(drg_img, 1727, 1386, 3784, 4345), TerminalProgressCallback("vw", "BLAH: "));
    
    write_image(orginal_tracks.c_str(), crop(drg_img, min_C, min_R, (max_C-min_C), (max_R - min_R)), TerminalProgressCallback("vw", "BLAH: "));
    printf("# cols = %d, # rows = %d\n",(max_C - min_C),(max_R - min_R));
        
    min_R =100000000;
    max_R = 0; 
    min_C = 100000000;
    max_C = 0;
    
    float holdX = 0.0;
    float holdY = 0.0;
    
    tracknum = 0;
    numpix = 0;
    point_size = 5;
    //typedef vector<vector<LOLAShot> >::const_iterator track_itr_t;
    //typedef vector<LOLAShot>::const_iterator shot_itr_t;
    //typedef vector<pointCloud>::const_iterator point_itr_t;
    for (track_itr_t track = track_list.begin(); track != track_list.end(); track++) {
      for (shot_itr_t shot = track->begin(); shot != track->end(); shot++) {
        for (point_itr_t point = shot->LOLAPt.begin(); point != shot->LOLAPt.end(); point++) {
          Vector2 pix = drg_georef.lonlat_to_pixel(SubVector<const Vector3>(point->coords, 0, 2));
          
          //convert pix location using affine transform assume:( x = pix(0), y = pix(1)) 
          holdX = pix.x();
          holdY = pix.y();
    
          pix.x() = d0*holdX + d1*holdY + d2;
          pix.y() = d3*holdX + d4*holdY + d5;
    
          
          if (BBox2(point_size, point_size, drg_img.cols() - point_size, drg_img.rows() - point_size).contains(pix)) {
            fill(crop(drg_img, int32(pix.x()) - point_size, int32(pix.y()) - point_size, point_size, point_size), PixelRGB<uint8>(255, 0, 0));
  
            //update current values for cropping
            row = int32( pix.y());
            col = int32( pix.x());
            Min_Max_Pixel( point_size, row, col, min_R, max_R, min_C, max_C);
 
            numpix++;
          }
        }
      }
      cout << "Status: " << ++tracknum << " / " << track_list.size() << "\r";
    }
     
    cout << endl << "numpix = " << numpix << endl;
    /*
    for (out_itr track; 
      for (int j = 0;

    } */
   
    //decide crop automatically


    write_image("hand_control_alignment.tif", crop(drg_img, min_C, min_R, (max_C-min_C), (max_R - min_R)), TerminalProgressCallback("vw", "BLAH moved tracks : "));
    //write_image(orginal_and_moved_tracks, crop(drg_img, min_C, min_R, (max_C-min_C), (max_R - min_R)), TerminalProgressCallback("vw", "BLAH moved tracks : "));
    printf("for tracks & aligned tracks # cols = %d, # rows = %d\n",(max_C - min_C),(max_R - min_R));

// write_image("final_really_big.tif",   crop(drg_img, 0, 0, 20000, 20000  )  , TerminalProgressCallback("vw", "big_save: "));
    return 0;

}


/*
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
*/















