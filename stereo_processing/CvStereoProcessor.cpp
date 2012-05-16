#include "CvStereoProcessor.h"
//#include "StereoParameters.h"

//#include "knShare/Image.h"

#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

//boost
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>


using namespace std;
using namespace cv;

//reads the calibration file
int readCalibrationFile(const string &calibrationFilename, 
                        cv::Mat &camera_matrix_left, cv::Mat &camera_matrix_right,
                        cv::Mat &dist_coeffs_left, cv::Mat &dist_coeffs_right, 
                        cv::Mat &R, cv::Mat &T, cv::Size &imageSize)
{
  ifstream calibrationFile (calibrationFilename.c_str());
  std::string line;
  double a00, a01,a02, a10, a11, a12, a20, a21, a22;
  int width, height;

  if (calibrationFile.is_open()){     
    std::string identifier;

    stringstream sline;  
    getline (calibrationFile,line);
    sline<<line;
    sline >> identifier >> a00 >> a01 >> a02 >> a10 >> a11>> a12 >> a20 >> a21 >> a22;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a20<<" "<<a21<<" "<<a22<<endl;
    camera_matrix_left.create(3, 3, CV_64F);
    camera_matrix_left.at<double>(0,0)=a00;
    camera_matrix_left.at<double>(0,1)=a01;
    camera_matrix_left.at<double>(0,2)=a02;
    camera_matrix_left.at<double>(1,0)=a10;
    camera_matrix_left.at<double>(1,1)=a11;
    camera_matrix_left.at<double>(1,2)=a12;
    camera_matrix_left.at<double>(2,0)=a20;
    camera_matrix_left.at<double>(2,1)=a21;
    camera_matrix_left.at<double>(2,2)=a22;
 
    stringstream sline1; 
    getline (calibrationFile,line);
    sline1 << line;
    sline1 >> identifier >> a00 >> a01 >> a02 >> a10;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<" "<<a10<<endl;
    dist_coeffs_left.create(4,1, CV_64F);
    dist_coeffs_left.at<double>(0,0)=a00;
    dist_coeffs_left.at<double>(1,0)=a01;
    dist_coeffs_left.at<double>(2,0)=a02;
    dist_coeffs_left.at<double>(3,0)=a10;

    stringstream sline2; 
    getline (calibrationFile,line);
    sline2 << line;
    sline2 >> identifier >> a00 >> a01 >> a02 >> a10 >> a11>> a12 >> a20 >> a21 >> a22;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a20<<" "<<a21<<" "<<a22<<endl;
    camera_matrix_right.create(3, 3, CV_64F);
    camera_matrix_right.at<double>(0,0)=a00;
    camera_matrix_right.at<double>(0,1)=a01;
    camera_matrix_right.at<double>(0,2)=a02;
    camera_matrix_right.at<double>(1,0)=a10;
    camera_matrix_right.at<double>(1,1)=a11;
    camera_matrix_right.at<double>(1,2)=a12;
    camera_matrix_right.at<double>(2,0)=a20;
    camera_matrix_right.at<double>(2,1)=a21;
    camera_matrix_right.at<double>(2,2)=a22;
  
    stringstream sline3; 
    getline (calibrationFile,line);
    sline3 << line;
    sline3 >> identifier >> a00 >> a01 >> a02 >> a10;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<" "<<a10<<endl;    
    dist_coeffs_right.create(4,1, CV_64F);
    dist_coeffs_right.at<double>(0,0)=a00;
    dist_coeffs_right.at<double>(1,0)=a01;
    dist_coeffs_right.at<double>(2,0)=a02;
    dist_coeffs_right.at<double>(3,0)=a10;

    stringstream sline4; 
    getline (calibrationFile,line);
    sline4 << line;
    sline4 >> identifier >> a00 >> a01 >> a02 >> a10 >> a11>> a12 >> a20 >> a21 >> a22;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<" "<<a10<<" "<<a11<<" "<<a12<<" "<<a20<<" "<<a21<<" "<<a22<<endl;
    R.create(3, 3, CV_64F);
    R.at<double>(0,0)=a00;
    R.at<double>(0,1)=a01;
    R.at<double>(0,2)=a02;
    R.at<double>(1,0)=a10;
    R.at<double>(1,1)=a11;
    R.at<double>(1,2)=a12;
    R.at<double>(2,0)=a20;
    R.at<double>(2,1)=a21;
    R.at<double>(2,2)=a22;


    stringstream sline5;    
    getline (calibrationFile,line);
    sline5 << line;
    sline5 >> identifier >> a00 >> a01 >> a02;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<endl;  
    T.create(3, 1, CV_64F);
    T.at<double>(0,0)=a00;
    T.at<double>(1,0)=a01;
    T.at<double>(2,0)=a02;

    stringstream sline6;
    getline (calibrationFile,line);
    sline6 << line;
    sline6 >> identifier >> width >> height;
    cout<<identifier<<" "<<width<<" "<<height<<endl;
    imageSize.width = width;
    imageSize.height = height;

    calibrationFile.close();
    return 1;
  }
  else {
    cout << "Unable to open file"; 
    return 0;
  }
}



CvStereoProcessor::CvStereoProcessor(CvStereoProcessorParameters const& params) :
    m_params(params),
    m_bMState(cvCreateStereoBMState())/*,
    m_imageDisparity(cvCreateMat(imageSize.height, imageSize.width, CV_16S)),
    m_imageDisparityNormalized(cvCreateMat(imageSize.height, imageSize.width, CV_8U)*/
  { 
    m_bMState->preFilterSize       = m_params.preFilterSize;
    m_bMState->preFilterCap        = m_params.preFilterCap;
    m_bMState->SADWindowSize       = m_params.sadWindowSize;
    m_bMState->minDisparity        = m_params.minDisparity;
    m_bMState->numberOfDisparities = m_params.numberOfDisparities;
    m_bMState->textureThreshold    = m_params.textureThreshold;
    m_bMState->uniquenessRatio     = m_params.uniquenessRatio;

    //rectification code - START
    //TODO:calibration filename can be a parameter in StereoCV21Parameters!!!
    string calibrationFilename="stereo_calibration.txt";
    cv::Mat camera_matrix[2], dist_coeffs[2];
    cv::Mat R, T;
    cv::Size imageSize;
    int readCalibFilenameError = readCalibrationFile(calibrationFilename,
                                                     camera_matrix[0], camera_matrix[1],
                                                     dist_coeffs[0], dist_coeffs[1],
                                                     R, T, imageSize);
    
    cv::Mat R1, R2, P1, P2;
    cv::Rect validRoi[2];

    cv::stereoRectify(camera_matrix[0], dist_coeffs[0],
                      camera_matrix[1], dist_coeffs[1],
                      imageSize, R, T,
                      R1, R2, P1, P2,
                      m_Q, 1,
                      cv::CALIB_ZERO_DISPARITY, imageSize, &validRoi[0], &validRoi[1]);

    cv::initUndistortRectifyMap(camera_matrix[0], dist_coeffs[0], R1, P1, imageSize, CV_16SC2, m_refRMap[0], m_refRMap[1]);
    cv::initUndistortRectifyMap(camera_matrix[1], dist_coeffs[1], R2, P2, imageSize, CV_16SC2, m_matchRMap[0], m_matchRMap[1]);
    //rectification code - END

    //the disparity map in OpenCV is stored on a 16bit format to increase the precision of the subpixel correlation
    //this effectively increases the disparity value by a factor of 16 and the statement below accopunts for this effect.
    m_Q.at<double>(3,2) = m_Q.at<double>(3,2)/16.0;
   m_Q.at<double>(3,2) =-m_Q.at<double>(3,2); //! @FIXME fix for OpenCV 2.3.1

   //TODO: copy m_Q to m_QOrig
   m_QOrig = m_Q; 
 
    m_imageDisparity=cvCreateMat(imageSize.height, imageSize.width, CV_16S);
    m_imageDisparityNormalized=cvCreateMat(imageSize.height, imageSize.width, CV_8U);
  }

#if 0
  CvStereoProcessor::CvStereoProcessor(CvStereoProcessorParameters const& params,
                                       cv::Mat const refRMap[2], 
                                       cv::Mat const matchRMap[2], 
                                       cv::Mat const& Q,
                                       cv::Size const& imageSize) :
    m_params(params),
    m_bMState(cvCreateStereoBMState()),
    //    m_refRMap(refRMap),
    //   m_matchRMap(matchRMap),
    m_Q(Q),
    m_imageDisparity(cvCreateMat(imageSize.height, imageSize.width, CV_16S)),
    m_imageDisparityNormalized(cvCreateMat(imageSize.height, imageSize.width, CV_8U))
  {
    // ???
    m_Q.at<double>(3,2) = m_Q.at<double>(3,2)/16.0;
    m_Q.at<double>(3,2) =-m_Q.at<double>(3,2); //fix for OpenCV 2.3.1

    
    //cv::Mat extRotationMat=cv:Matt::eye(4,4, CV_64F);
    //for(i = 0; i < 3; i++){
    //  for(j =0; j < 3; j++){
    //     extRotationMat.at<double>(i,j) = rotationMat.at<double>(i,j);
    //	}
    //}	
    //TODO: multiply m_Q = extRotationMatrix*m_Q;
   

    // ???
    // dunno why this does not work as member-ctor

    //the rectification maps for the reference and match images
    m_refRMap[0] = refRMap[0];
    m_refRMap[1] = refRMap[1];
    m_matchRMap[0] = matchRMap[0];
    m_matchRMap[1] = matchRMap[1];

    m_bMState->preFilterSize       = m_params.preFilterSize;
    m_bMState->preFilterCap        = m_params.preFilterCap;
    m_bMState->SADWindowSize       = m_params.sadWindowSize;
    m_bMState->minDisparity        = m_params.minDisparity;
    m_bMState->numberOfDisparities = m_params.numberOfDisparities;
    m_bMState->textureThreshold    = m_params.textureThreshold;
    m_bMState->uniquenessRatio     = m_params.uniquenessRatio;

  }
#endif
  CvStereoProcessor::~CvStereoProcessor()
  {
    cvReleaseMat(&m_imageDisparity);
    cvReleaseMat(&m_imageDisparityNormalized);
    cvReleaseStereoBMState(&m_bMState);
  }

  void
  CvStereoProcessor::savePoints(string const &filenameNoPath)
  {
    double const maxZ = 1.0e4;
     
    ofstream ostr;
    string pointCloudFilename = m_params.resDir+filenameNoPath;
    ostr.open(pointCloudFilename.c_str(), ios::out);    

    for (int y = 0; y < m_points.rows; ++y) {
      for (int x = 0; x < m_points.cols; ++x) {
        
        Vec3f const& point = m_points.at<Vec3f>(y, x);
        //if (!(fabs(point[2] - maxZ) < FLT_EPSILON || fabs(point[2]) > maxZ)) {
        if (fabs(point[2]-maxZ)> FLT_EPSILON && fabs(point[2]<maxZ) && point[2]>0){
          ostr << point[0] << " " << point[1] << " " << point[2] << endl;
        }
      }
    }
    ostr.close();
  }

#if 0
  void 
  CvStereoProcessor::iplImageFromKnImage(IplImage& ipl, Image const& img)
  {
    ipl.ID = 0;              /* version (=0)*/
    ipl.nChannels = 1;      /* Most of OpenCV functions support 1,2,3 or 4 channels */
    //    ipl.alphaChannel;      /* Ignored by OpenCV */

    // this could be more generic, but that's a case-statement...
    ipl.depth = IPL_DEPTH_8U; /* Pixel depth in bits: IPL_DEPTH_8U, IPL_DEPTH_8S, IPL_DEPTH_16S,
                               IPL_DEPTH_32S, IPL_DEPTH_32F and IPL_DEPTH_64F are supported.  */
    //    char colorModel[4];     /* Ignored by OpenCV */
    //    char channelSeq[4];     /* ditto */
    ipl.dataOrder = 0;       /* 0 - interleaved color channels, 1 - separate color channels.
                                cvCreateImage can only create interleaved images */
    ipl.origin= 0;          /* 0 - top-left origin,
                               1 - bottom-left origin (Windows bitmaps style).  */
    //    int  align;             /* Alignment of image rows (4 or 8).
    //                               OpenCV ignores it and uses widthStep instead.    */
    ipl.width = img.width(); /* Image width in pixels.                           */
    ipl.height = img.height();            /* Image height in pixels.                          */
    ipl.roi = NULL;              /* Image ROI. If NULL, the whole image is selected. */
    ipl.maskROI = NULL;      /* Must be NULL. */
    ipl.imageId = NULL;               /* "           " */
    ipl.tileInfo = NULL;  /* "           " */
    ipl.imageSize = ipl.width * ipl.height;         /* Image data size in bytes
                               (==image->height*image->widthStep
                               in case of interleaved data)*/
    ipl.imageData = (char *)img.data();        /* Pointer to aligned image data.         */
    ipl.widthStep = ipl.width;         /* Size of aligned image row in bytes.    */
    //    int  BorderMode[4];     /* Ignored by OpenCV.                     */
    //    int  BorderConst[4];    /* Ditto.                                 */
    ipl.imageDataOrigin = NULL;  /* Pointer to very origin of image data
                               (not necessarily aligned) -
                               needed for correct deallocation */
  }
#endif

/*
  //this function should be used but for now is not working
  CvMat const&
  CvStereoProcessor::normalizedDisparityImage()
  {
    cvNormalize(m_imageDisparity, m_imageDisparityNormalized, 0, 256, CV_MINMAX);
    return *m_imageDisparityNormalized;
  }
*/


  //modifies the Q matrix with the rotation  and translation
  void CvStereoProcessor::changeCoordinates(cv::Mat const &rotationMat, cv::Mat const &vectorMat)
  {
    //define the extended rotation matrix
    cv::Mat extRotationMat=cv::Mat::eye(4,4, CV_64F);
    for(int i = 0; i < 3; i++){
      for(int j =0; j < 3; j++){
         extRotationMat.at<double>(i,j) = rotationMat.at<double>(i,j);
      }
    } 

    for (int i = 0; i < 3; ++i) {
      extRotationMat.at<double>(i,3) = vectorMat.at<double>(i);
    }

    //multiply m_Q = extRotationMatrix*m_Q;
    int flags = 0;    
    gemm(extRotationMat, m_QOrig, 1.0, /*NULL*/ m_Q, 0.0, m_Q, flags);
   
  }

  void
  CvStereoProcessor::processImagePair(IplImage *refImg, IplImage *matchImg)
  {
    Mat refMat = refImg;
    Mat matchMat = matchImg;

    //image rectification before stereo correlation
    remap(refMat, m_refImgRect, m_refRMap[0], m_refRMap[1], CV_INTER_LINEAR);
    remap(matchMat, m_matchImgRect, m_matchRMap[0], m_matchRMap[1], CV_INTER_LINEAR);
    
    IplImage refImgRect;
    IplImage matchImgRect;
    refImgRect = m_refImgRect;
    matchImgRect = m_matchImgRect;

    cvFindStereoCorrespondenceBM(&refImgRect, &matchImgRect, m_imageDisparity, m_bMState);

    // Ara's way to do the conversion
    //Mat dispMat = &(*imageDisparity);

    Mat disMat = m_imageDisparity;

    reprojectImageTo3D(disMat, m_points, m_Q, 1);
  }
/*
  void
  CvStereoProcessor::processImagePair(kn::Image const& refImg, kn::Image const& matchImg)
  {
    IplImage ref;
    IplImage match;

    iplImageFromKnImage(ref, refImg);
    iplImageFromKnImage(match, matchImg);

    processImagePair(&ref, &match);
  }
*/
  //save the normalized disparity for display and debug  
  void
  CvStereoProcessor::saveDisparity(string const& filenameNoPath)
  {
    string disparityFilename = m_params.resDir+filenameNoPath;
    cvNormalize(m_imageDisparity, m_imageDisparityNormalized, 0, 256, CV_MINMAX);
    cvSaveImage( disparityFilename.c_str(), m_imageDisparityNormalized);
  }
  


