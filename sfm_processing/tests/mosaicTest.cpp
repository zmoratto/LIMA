#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

// Program includes
#include "../mosaic.h"

int main(int argc, char* argv[])
{
  string imgDir, pcDir, resultsDir, configFilename;
  
  if(argc!=4){
     cout<<"Usage: ./mosaic imgDir pcDir resultsDir"<<endl;
     cout<<"Requires list of stereo pairs in dataDir named test-stereo-pairs.txt"<<endl;
     cout<<"Using default: ./mosaic mosaicTestData mosaicResults"<<endl;
     sleep(5);
     imgDir = string("mosaicTestData") +string("/");
     resultsDir = string("mosaicResults") +string("/");
  }
  else {
     imgDir = string(argv[1])+string("/");
     pcDir = string(argv[2])+string("/");
     resultsDir = string(argv[3]) + string("/");
  }

  configFilename = "mosaic_settings.txt";
  mosaicSettings settings = ReadMosaicSettings(configFilename);
  mosaicProcessor mosaic(settings);
  PrintSettings(settings);
  cout<<"imgDir="<<imgDir<<endl;
  cout<<"pcDir="<<pcDir<<endl;
  cout<<"resultsDir="<<resultsDir<<endl;

  //create the results directory
  string command = string("mkdir ") + resultsDir;
  cout<<command<<endl;
  system(command.c_str());  

  int numPairs;
  vector<ImagePairNames> imgPairVector;
  //vector<CameraPairNames> camPairVector;

  string leftCameraModelName, rightCameraModelName;
  
  string listFilename = imgDir + string("test-stereo-pairs.txt");
  
  // Try to open the image list file
  ifstream file;
  file.open(listFilename.c_str());
  if (file.is_open())
  {
     file.close();
     imgPairVector = mosaic.ReadListFile(listFilename);
  }
  else
  {
    cout<<"list file "<<listFilename<<" not found. exiting."<<endl;
    return(0);

  }

  BBox mosaicBBox;
  vector<BBox> BBoxArray;
  float radPerPixX, radPerPixY;

  //calculates the boundaries of the mosaic
  cout<<"computing the boundary boxes"<<endl;
  mosaic.calcBoundaryBoxes( imgDir, imgPairVector, BBoxArray, mosaicBBox, radPerPixX, radPerPixY );

  cout<<"mosaic bounding box:"<< "minTilt= "<<mosaicBBox.minTilt<<", maxTilt= "<< mosaicBBox.maxTilt<<
                                ", minPan= "<< mosaicBBox.minPan<<", maxPan= "<<mosaicBBox.minPan<<endl;

  //compute the image mosaic without tiling. this function will be removed. for now it is only for validation of the results
  if (settings.makeImageMosaic == 1){
    cout<<"display corrected image mosaic"<<endl;
    mosaic.displayCorrectedImageMosaic(resultsDir, imgDir, imgPairVector, BBoxArray, mosaicBBox, radPerPixX, radPerPixY );
  }
  else{
    cout<<"skipping display corrected mosaic"<<endl;
  }

  //compute the tiles of the point cloud and image mosaic
  cout<<"creating image and point cloud tiles ..."<<endl;
  mosaic.makeTiles(imgDir, pcDir, resultsDir, imgPairVector, BBoxArray, mosaicBBox, radPerPixX, radPerPixY);

  //compute the composite mosaic from all existing tiles
  if (settings.makeTileMosaic == 1){
    cout<<"display tile mosaic from tiles"<<endl;
    mosaic.displayTileMosaic(resultsDir,mosaicBBox, radPerPixX, radPerPixY);
    cout<<"display image mosaic from tiles"<<endl;
    mosaic.displayPCMosaic(resultsDir,mosaicBBox, radPerPixX, radPerPixY);
  }

  return 0;
}
