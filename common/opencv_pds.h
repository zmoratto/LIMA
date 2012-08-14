#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "pds_read.h"

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

#include <sys/time.h>

using namespace std;
using namespace cv;

int readPDSFileToIplImage( const string &filename, IplImage* &dest );
int readPDSFileToMat( const string &filename, Mat &dest );

class PDSTileReader
   {
   public:
   PDSTileReader();
   ~PDSTileReader();
   int initialize(const string &inputFilename);
   int getHeaderInfo( long &cols, long &rows, long &bands, long &bits);
   int readPDSTileToMat( Mat &dest, CvRect roi);
   int readPDSTileToIplImage( IplImage* &dest, CvRect roi);

   private:
   unsigned short maxVal;
   long nRows;
   long nCols;
   long nBands;
   long nBits;
   ODLTREE odltree;
   ODLTREE image_node;
   bool initialized;
   string filename;
   int _readPartialPDSImage( unsigned char* &dataBuffer, CvRect roi );
};

int _readFullPDSImage(const string &strInFilename, unsigned char* &dataBuffer, 
                   long &nCols, long &nRows, long &nBands, long &nBits);

