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
  string dataDir, resultsDir, configFilename;
  if(argc!=3)
     {
     cout<<"Usage: ./mosaic dataDir resultsDir"<<endl;
     cout<<"Requires list of stereo pairs in dataDir named test-stereo-pairs.txt"<<endl;
     }

  dataDir = string(argv[1])+string("/");
  resultsDir = string(argv[2]) + string("/");
  configFilename = "mosaic_settings.txt";
  mosaicSettings settings = ReadMosaicSettings(configFilename);
  mosaicProcessor mosaic(settings);
  PrintSettings(settings);
  cout<<"dataDir="<<dataDir<<endl;
  cout<<"resultsDir="<<resultsDir<<endl;

  //create the results directory
  string command = string("mkdir ") + resultsDir;
  cout<<command<<endl;
  system(command.c_str());  

  int numPairs;
  vector<ImagePairNames> imgPairVector;
  //vector<CameraPairNames> camPairVector;

  string leftCameraModelName, rightCameraModelName;
  
  string listFilename = dataDir + string("test-stereo-pairs.txt");
  
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

  mosaic.calcBoundaryBoxes( dataDir, imgPairVector, BBoxArray, mosaicBBox, radPerPixX, radPerPixY );

  if (settings.makeImageMosaic == 1){
    mosaic.displayCorrectedImageMosaic(resultsDir, dataDir, imgPairVector, BBoxArray, mosaicBBox, radPerPixX, radPerPixY );
  }

  mosaic.makeTiles(dataDir, resultsDir, imgPairVector, BBoxArray, mosaicBBox, radPerPixX, radPerPixY);

  if (settings.makeTileMosaic == 1){
    mosaic.displayTileMosaic(resultsDir);
  }

  return 0;
}
