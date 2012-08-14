#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "pds_read.h"
#include "opencv_pds.h"

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

#include <sys/time.h>

using namespace std;
using namespace cv;


int readPDSFileToIplImage( const string &filename, IplImage* &dest )
{

uchar* data;
long nCols, nRows, nBands, nBits;

//read pds image
if( _readFullPDSImage(filename, data, nCols, nRows, nBands, nBits) ){
   cout<<"Error: image not loaded"<<endl;
   return -1;
   }

//check number of bands
if(nBands != 1 && nBands !=3){
   cout<<"Error: only nBands=1 or nBands=3 supported" << endl;
   return -1;
   }

// create image
dest = cvCreateImage(Size(nCols, nRows), IPL_DEPTH_8U, nBands);

// get pointer to image data
uchar* ptr = (uchar*) dest->imageDataOrigin;
memcpy(ptr, data, nRows*nCols*nBands);

// release data buffer
delete[] data;

 return 0;
}

int readPDSFileToMat( const string &filename, Mat &dest )
{

uchar* data;
long nCols, nRows, nBands, nBits;

//read pds image
if( _readFullPDSImage(filename, data, nCols, nRows, nBands, nBits) ){
   cout<<"Error: image not loaded"<<endl;
   return -1;
   }

//check number of bands
if(nBands != 1 && nBands !=3){
   cout<<"Error: only nBands=1 or nBands=3 supported" << endl;
   return -1;
   }

// create image
dest.create(nRows, nCols, CV_MAKETYPE(CV_8U, nBands) );

// get pointer to image data
uchar* ptr = (uchar*) dest.ptr(0);
memcpy(ptr, data, nRows*nCols*nBands);

// release data buffer
delete[] data;
   return 0;
}

// allocates memory and stores image in dataBuffer
int _readFullPDSImage(const string &strInFilename, unsigned char* &dataBuffer, 
                   long &nCols, long &nRows, long &nBands, long &nBits)
{

    // parse the pds file to get the top level node odltree
   ODLTREE odltree = OaParseLabelFile((char *)strInFilename.c_str(), 
				       (char*) ERRS_LOC, ODL_EXPAND_STRUCTURE, TRUE);

    if (odltree == NULL) { 
        OaReportError( (char*) "Error in OaParseLabelFile()" );
        return -1; 
    }

    // grab the image node from the top level node
    ODLTREE image_node = OdlFindObjDesc(odltree, (char*) "IMAGE", NULL, NULL,
                                        ODL_RECURSIVE_DOWN, ODL_TO_END);

    if (image_node == NULL) {
        OaReportError( (char*) "Error in OaParseLabelFile()" );
        return -1;
    }

    ///////////////////////////////////////////////////////////////////////////
    // process label part of PDS file
    ///////////////////////////////////////////////////////////////////////////

    nBands = get_pds_label_value(image_node, "BANDS");
    nBits =  get_pds_label_value(image_node, "SAMPLE_BITS");
    if(nBits != 8 && nBits != 16 && nBits != 32){
        OaReportError( (char*) "Error only 8, 16, and 32 bit images supported" );
        return -1;
        }

    nRows = get_pds_label_value(image_node, "LINES");
    nCols = get_pds_label_value(image_node, "LINE_SAMPLES");

    ///////////////////////////////////////////////////////////////////////////
    // if data is 16 bit we need the max value to convert to 8 bit
    ///////////////////////////////////////////////////////////////////////////
    double factor = 1.0;
    if(nBits == 16)
       {
       unsigned short maxVal=0;
       for(int band=0; band<nBands; ++band)
          {
          // open the image (the bands are 1 indexed here)
          OA_OBJECT image_handle = OaOpenImage(image_node, band+1);
          if (image_handle == NULL) {
              OaReportError( (char*) "Error in OaOpenImage()" );
              return -1;
              }
          OA_OBJECT image_obj = OaReadImage(image_node, band+1);
          if (image_obj == NULL) {
              OaReportError( (char*) "Error in OaReadImage()" );
              return -1;
              }

          unsigned short* sptr = (unsigned short*) image_obj->data_ptr;
          maxVal = max(*max_element(sptr, sptr+nRows*nCols), maxVal);

          // free image obj memory
          OaCloseImage(image_handle);
          OaDeleteObject(image_obj);
          }

       if(maxVal!=0)
          factor = 256.0/maxVal;
       }


    ///////////////////////////////////////////////////////////////////////////
    // process data part of PDS file (image data)
    ///////////////////////////////////////////////////////////////////////////

    // allocate memory
    dataBuffer = new unsigned char[nRows*nCols*nBands];

    for(int band=0; band<nBands; ++band)
         {
         int opencvBand;
         if(nBands == 3)
            opencvBand = 2-band;
         else 
            opencvBand = band;

         // open the image (the bands are 1 indexed here)
         OA_OBJECT image_handle = OaOpenImage(image_node, band+1);
         if (image_handle == NULL) {
             OaReportError( (char*) "Error in OaOpenImage()" );
             return -1;
             }

         // NOTE: batch version -- no copying done, you may want to read the
         // image by chunks here, if it's huge, using OaReadPartialImage()
         OA_OBJECT image_obj = OaReadImage(image_node, band+1);
         if (image_obj == NULL) {
             OaReportError( (char*) "Error in OaReadImage()" );
             return -1;
             }

         if(nBits == 8)
               {
               unsigned char* uptr = (unsigned char*)image_obj->data_ptr;
               for(long i = 0; i < nRows * nCols; ++i)
                  dataBuffer[i*nBands+opencvBand] = uptr[i];
               }
         else if(nBits == 16)
               {
               short unsigned int* sptr = (short unsigned int*) image_obj->data_ptr;
               for(long i = 0; i < nRows * nCols; ++i)
                  dataBuffer[i*nBands+opencvBand] = (uchar) (sptr[i] * factor);          
               }
         else if(nBits == 32)
               {
               float* fptr = (float*) image_obj->data_ptr;
               for(long i = 0; i < nRows * nCols; ++i)
                  dataBuffer[i*nBands+opencvBand] = (uchar) (fptr[i] * 256.0);          
               }

         // free image obj memory
         OaCloseImage(image_handle);
         OaDeleteObject(image_obj);

         }


    // free tree memory
    OdlFreeTree( odltree );

    return 0;
}

PDSTileReader::PDSTileReader()
{
   initialized = false;
}

PDSTileReader::~PDSTileReader()
{
   if(initialized)
      OdlFreeTree( odltree );
}

int PDSTileReader::initialize(const string &inputFilename)
{
   if(initialized)
      OdlFreeTree( odltree );

   filename = inputFilename;
   maxVal = 0;
   // parse the pds file to get the top level node odltree
   odltree = OaParseLabelFile((char *)filename.c_str(), 
				   (char*) ERRS_LOC, ODL_EXPAND_STRUCTURE, TRUE);

   if (odltree == NULL) { 
      OaReportError( (char*) "Error in OaParseLabelFile()" );
      return -1; 
   }

   // grab the image node from the top level node
   image_node = OdlFindObjDesc(odltree, (char*) "IMAGE", NULL, NULL,
                               ODL_RECURSIVE_DOWN, ODL_TO_END);

   if (image_node == NULL) {
      OaReportError( (char*) "Error in OaParseLabelFile()" );
      return -1;
   }


   ///////////////////////////////////////////////////////////////////////////
   // process label part of PDS file
   ///////////////////////////////////////////////////////////////////////////

   nBands = get_pds_label_value(image_node, "BANDS");
   nBits =  get_pds_label_value(image_node, "SAMPLE_BITS");
   if(nBits != 8 && nBits != 16 && nBits != 32){
       OaReportError( (char*) "Error only 8, 16, and 32 bit images supported" );
       return -1;
       }

   nRows = get_pds_label_value(image_node, "LINES");
   nCols = get_pds_label_value(image_node, "LINE_SAMPLES");


   ///////////////////////////////////////////////////////////////////////////
   // process data part of PDS file (image data)
   ///////////////////////////////////////////////////////////////////////////

   // if the data is 16 bit we need the max value
   if(nBits == 16)
      {

      short int * buffer = new short int[nCols];

      for(int band=0; band<nBands; ++band)
         {
         // open the image (the bands are 1 indexed here)
         OA_OBJECT image_handle = OaOpenImage(image_node, band+1);
         if (image_handle == NULL) {
             OaReportError( (char*) "Error in OaOpenImage()" );
             return -1;
             }
         
         // read image by rows
         for(int row = 0; row < nRows; ++row)
            {
            // read image (rows and cols are 1 indexed)
            OA_OBJECT image_obj = OaReadPartialImage(image_handle, row+1, row+1,
                                               1, nCols);
            if (image_obj == NULL) {
                OaReportError( (char*) "Error in OaReadImage()" );
                return -1;
                }

            unsigned short* ptr = (unsigned short*) image_obj->data_ptr;
            
            // get new max
            maxVal = max(*max_element(ptr, ptr+nCols), maxVal);

            // delete object
            OaDeleteObject(image_obj);
            }

         // close image
         OaCloseImage(image_handle);   

         }
      
      // release line buffer
      delete[] buffer;

      }

   // set initialized to true
   initialized = true;

   return 0;
}


int PDSTileReader::getHeaderInfo( long &cols, long &rows, long &bands, long &bits)
{
   // check if initialized
   if( !initialized ){
      cout<<"Error: PDSTileReader not initialized"<<endl;
      return -1;
      }  

   cols = nCols;
   rows = nRows;
   bands = nBands;
   bits = nBits;
   return 0;
}

int PDSTileReader::readPDSTileToIplImage( IplImage* &dest, CvRect roi)
{

uchar* data;

// check if initialized
if( !initialized ){
   cout<<"Error: PDSTileReader not initialized"<<endl;
   return -1;
   }  

//check number of bands
if(nBands != 1 && nBands !=3){
   cout<<"Error: only nBands=1 or nBands=3 supported" << endl;
   return -1;
   }

//check tile parameters
if( roi.x+roi.width > nCols || roi.x < 0 || roi.y+roi.height > nRows 
    || roi.y < 0 || roi.width <= 0 || roi.height <= 0 ){
   cout<<"Error: tile parameters not valid" << endl;
   return -1;
   }

if(_readPartialPDSImage(data, roi) ){
   cout<<"Error: image not loaded"<<endl;
   return -1;
   }  

// create image
dest = cvCreateImage(Size(roi.width, roi.height), IPL_DEPTH_8U, nBands);
// get pointer to image data
uchar* ptr = (uchar*) dest->imageDataOrigin;
// copy data to image
memcpy(ptr,data,roi.width*roi.height*nBands);
// release data buffer
delete[] data;

return 0;

}


int PDSTileReader::readPDSTileToMat( Mat &dest, CvRect roi)
{
uchar* data;

// check if initialized
if( !initialized ){
   cout<<"Error: PDSTileReader not initialized"<<endl;
   return -1;
   }  

//check number of bands
if(nBands != 1 && nBands !=3){
   cout<<"Error: only nBands=1 or nBands=3 supported" << endl;
   return -1;
   }

//check tile parameters
if( roi.x+roi.width > nCols || roi.x < 0 || roi.y+roi.height > nRows 
    || roi.y < 0 || roi.width <= 0 || roi.height <= 0 ){
   cout<<"Error: tile parameters not valid" << endl;
   return -1;
   }

if(_readPartialPDSImage(data, roi) ){
   cout<<"Error: image not loaded"<<endl;
   return -1;
   }  

// create image
dest.create(roi.height, roi.width, CV_MAKETYPE(CV_8U, nBands) );
// get pointer to image data
uchar* ptr = (uchar*) dest.ptr(0);
// copy data to image
memcpy(ptr,data,roi.width*roi.height*nBands);
// release data buffer
delete[] data;

return 0;
}


int PDSTileReader::_readPartialPDSImage( unsigned char* &dataBuffer, CvRect roi )
{
   // check if initialized 
   if( !initialized ){
      cout<<"Error: PDSTileReader not initialized"<<endl;
      return -1;
      }  

    if(nBits != 8 && nBits != 16){
        OaReportError( (char*) "Error only 8 bit and 16 bit images supported" );
        return -1;
        }

    ///////////////////////////////////////////////////////////////////////////
    // process data part of PDS file (image data)
    ///////////////////////////////////////////////////////////////////////////

    // allocate memory
    dataBuffer = new unsigned char[roi.height*roi.width*nBands];

    for(int band=0; band<nBands; ++band)
         {
         int opencvBand;
         if(nBands == 3)
            opencvBand = 2-band;
         else 
            opencvBand = band;

         // open the image (the bands are 1 indexed here)
         OA_OBJECT image_handle = OaOpenImage(image_node, band+1);
         if (image_handle == NULL) {
             OaReportError( (char*) "Error in OaOpenImage()" );
             return -1;
             }

         // read image tile (the +1 is to account for 1 indexing)
         OA_OBJECT image_obj = OaReadPartialImage(image_handle, roi.y+1, roi.y+roi.height,
                                               roi.x+1, roi.x+roi.width);
         if (image_obj == NULL) {
             OaReportError( (char*) "Error in OaReadImage()" );
             return -1;
             }

         if(nBits == 8)
               {
               unsigned char* uptr = (unsigned char*)image_obj->data_ptr;
               for(long i = 0; i < roi.height * roi.width; ++i)
                  dataBuffer[i*nBands+opencvBand] = uptr[i];
               }
         else if(nBits == 16)
               {
               double factor = 1.0;
               if( maxVal>0 )
                  factor = 256.0/maxVal;
               short unsigned int* sptr = (short unsigned int*) image_obj->data_ptr;
               for(long i = 0; i < roi.height * roi.width; ++i)
                  dataBuffer[i*nBands+opencvBand] = (uchar) (sptr[i] * factor);          
               }
         else if(nBits == 32)
               {
               float* fptr = (float*) image_obj->data_ptr;
               for(long i = 0; i < roi.height * roi.width; ++i)
                  dataBuffer[i*nBands+opencvBand] = (uchar) (fptr[i] * 256.0);          
               }


         // free image obj memory
         OaCloseImage(image_handle);
         OaDeleteObject(image_obj);

         }

    return 0;
}

