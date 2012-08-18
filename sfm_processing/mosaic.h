// mosaic.cpp : mosaicing images and point cloud using tiles

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

using namespace std;

struct ImagePairNames{
std::string left;
std::string right;
};

struct point3D
{
  float x;
  float y;
  float z;
};

struct BBox
{
  float minPan;
  float minTilt;
  float maxPan;
  float maxTilt;
};

struct mosaicSettings
{
  double pixelSize;// = 0.000012;
  int tileWidth;// = 512;
  int tileHeight;// = 512;
  float scaleFactor;// = 0.25;
  int makeTileMosaic;// = 1;
  int makeImageMosaic;//= 1;
};


// function to read the settings
struct mosaicSettings ReadMosaicSettings(const string &settingsFilename);

void PrintSettings(struct mosaicSettings settings);

class mosaicProcessor
{
public:
   //constructor for passing settings
   mosaicProcessor(const mosaicSettings &params);

   // function to determine which images overlap the tile
   std::vector<int> makeOverlapListFromBB(struct BBox tileBox, std::vector<BBox> imageBox);

   // function to display the tile mosaic for debug
   void displayTileMosaic(string resultsDir, BBox mosaicBBox, double radPerPixX, double radPerPixY);

   //function to create the point cloud mosaic
   void displayPCMosaic(string resultsDir, BBox mosaicBBox, double radPerPixX, double radPerPixY);

   //reads a list file and creates a vector of image pairs
   //vector<struct imgPair> ReadListFile(string listFilename);
   //reads a list file and creates a vector of image pairs
   vector<ImagePairNames> ReadListFile(string listFilename);

   std::vector<point3D> readXYZFile(const string &inputFilename, int rows, int cols);

   void saveXYZFile(const string &inputFilename, std::vector<point3D> &points, unsigned char *uchar_tile);

   // function to read the cam params from a .cahv file
   void ReadCamParamsFromFile(char *filename, double*c, double*a, double*h, double*v, int* r_rows, int* r_cols);

   // function to take a quaternion and compute the 3x3 rotation matrix
   void quaternionToRotation( double quaternion[4], double matrix[3][3] );

   //rotates the input vector by the quaternion
   void rotateWithQuaternion( double input[3], double output[3], double quaternion[4] );

   // function to compute the dot product of 3d vectors
   inline double dotprod3( double A[3], double B[3] );

   // function to get the x and y positions in an image from the 3d point p
   void getXYFromCAHV( double c[3], double a[3], double h[3], double v[3], double &x, double &y, double p[3] );

   // function to read cam params from a PDS file
   void ReadCamParamsFromPDSFile(char *filename, double*c, double*a, double*h, double*v, int* r_rows, int* r_cols);

   // function to display the image mosaic corrected using the quaternion
   void displayCorrectedImageMosaic(string resultsDir, string dataDir, vector<ImagePairNames> imgPairVector, 
                                    vector<BBox> BBoxArray, struct BBox mosaicBBox,  float radPerPixX, 
                                    float radPerPixY);

   // original function to display the image mosaic
   void displayImageMosaic(string resultsDir, string dataDir, vector<ImagePairNames> imgPairVector, 
                           struct BBox mosaicBBox,  float radPerPixX, float radPerPixY);

   // function that makes the tiles
   void makeTiles(string dataDir, string resultsDir, std::vector<ImagePairNames> imgPairVector, std::vector<BBox> BBoxArray, 
                  struct BBox mosaicBBox, double radPerPixX, double radPerPixY);

   //function to calculate boundary boxes for each image pair
   void calcBoundaryBoxes( string dataDir, vector<ImagePairNames> &imgPairVector, vector<BBox> &BBoxArray, BBox &mosaicBBox, float &radPerPixX, float &radPerPixY );

protected: 
  double pixelSize;
  int tileWidth;
  int tileHeight;
  float scaleFactor;
};



