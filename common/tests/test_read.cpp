#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "opencv_pds.cpp"

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
Mat tile;
IplImage *tile2;
CvRect roi;
roi.x=0;roi.y=0;roi.width=600;roi.height=600;

// initialize tile reader
tileReader.init(filename);
long nCols, nRows, nBands, nBits;
tileReader.getHeaderInfo(nCols,nRows,nBands,nBits);
cout<<"bands: "<<nBands<<" bits: "<<nBits<<endl;

// read tile to IplImage
if(tileReader.readPDSTileToIplImage(tile2,roi))
   return -1;
cvSaveImage("tileIpl.pgm",tile2);
cvReleaseImage(&tile2);

// read tile to Mat
if(tileReader.readPDSTileToMat(tile,roi))
   return -1;
imwrite("tileMat.pgm",tile);



IplImage *im2;

// read file to IplImage
if(readPDSFileToIplImage(filename,im2))
   return -1;
cvSaveImage( "iplImage.pgm", im2);
cvReleaseImage(&im2);


Mat im;
// read file to Mat
if(readPDSFileToMat(filename,im))
   return -1;
imwrite("mat.pgm",im);


 return 0;
}



