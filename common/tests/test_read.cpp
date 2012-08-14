#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "../opencv_pds.h"

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

using namespace std;
using namespace cv;

int main(int argc, char **argv )
{

//check arg num
if(argc!=2){
   cout<<"Usage: \n"<<"pds_read <filename>"<<endl; 
   return 0;
   }

string filename = argv[1];

// declarations for tile reader
PDSTileReader tileReader;
IplImage *tile;

// initialize tile reader
tileReader.initialize(filename);
long nCols, nRows, nBands, nBits;
tileReader.getHeaderInfo(nCols,nRows,nBands,nBits);
cout<<"bands: "<<nBands<<" bits: "<<nBits<<endl;

// read tile to IplImage
CvRect roi;
roi.x=0;roi.y=0;roi.width=nCols/2;roi.height=nCols/2;
if(tileReader.readPDSTileToIplImage(tile,roi))
   return -1;
cvSaveImage("test_tile.pgm",tile);
cvReleaseImage(&tile);

// read file to IplImage
IplImage *im;
if(readPDSFileToIplImage(filename,im))
   return -1;
cvSaveImage( "iplImage.pgm", im);
cvReleaseImage(&im);

 return 0;
}



