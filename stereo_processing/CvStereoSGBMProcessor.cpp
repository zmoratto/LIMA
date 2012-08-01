#include "CvStereoSGBMProcessor.h"
//#include "StereoParameters.h"

//#include "knShare/Image.h"

#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

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



CvStereoSGBMProcessor::CvStereoSGBMProcessor(CvStereoSGBMProcessorParameters const& params) :
    m_params(params)/*,
    m_imageDisparity(cvCreateMat(imageSize.height, imageSize.width, CV_16S)),
    m_imageDisparityNormalized(cvCreateMat(imageSize.height, imageSize.width, CV_8U)*/
  { 
    sgbm.preFilterCap              = m_params.preFilterCap;
    sgbm.SADWindowSize             = m_params.sadWindowSize;
    sgbm.P1                        = m_params.P1;
    sgbm.P2                        = m_params.P2;
    sgbm.minDisparity              = m_params.minDisparity;
    sgbm.numberOfDisparities       = m_params.numberOfDisparities;
    sgbm.uniquenessRatio           = m_params.uniquenessRatio;
    NEED_RECTIFICATION             = m_params.needRectification;
    sgbm.fullDP=true;

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

    // get new size and cam matrix for new size
    smallSize.height = round(imageSize.height*m_params.scaleFactor);
    smallSize.width = round(imageSize.width*m_params.scaleFactor);
    camera_matrix[0]=camera_matrix[0]*m_params.scaleFactor;
    camera_matrix[1]=camera_matrix[1]*m_params.scaleFactor;


    if(NEED_RECTIFICATION)
    {
      //FIXED: order of cv::CALIB_ZERO_DISPARITY, 1 was switched
      cv::stereoRectify(camera_matrix[0], dist_coeffs[0],
                        camera_matrix[1], dist_coeffs[1],
                        smallSize, R, T,
                        R1, R2, P1, P2,
                        m_Q, cv::CALIB_ZERO_DISPARITY,
                        1, smallSize, &validRoi[0], &validRoi[1]);

      cv::initUndistortRectifyMap(camera_matrix[0], dist_coeffs[0], R1, P1, smallSize, CV_16SC2, m_refRMap[0], m_refRMap[1]);
      cv::initUndistortRectifyMap(camera_matrix[1], dist_coeffs[1], R2, P2, smallSize, CV_16SC2, m_matchRMap[0], m_matchRMap[1]);

      //the disparity map in OpenCV is stored on a 16bit format to increase the precision of the subpixel correlation
      //this effectively increases the disparity value by a factor of 16 and the statement below accopunts for this effect.
      m_Q.at<double>(3,2) = m_Q.at<double>(3,2)/16.0;
      m_Q.at<double>(3,2) =-m_Q.at<double>(3,2); //! @FIXME fix for OpenCV 2.3.1

      //TODO: copy m_Q to m_QOrig
      m_QOrig = m_Q; 
      //cout<<m_Q<<endl; */
      //rectification code - END
      }

    m_imageDisparity=cvCreateMat( smallSize.height, smallSize.width, CV_16S);
    m_imageDisparityNormalized=cvCreateMat( smallSize.height, smallSize.width, CV_8U);
  }

  CvStereoSGBMProcessor::~CvStereoSGBMProcessor()
  {
  }

  void
  CvStereoSGBMProcessor::savePoints(string const &filenameNoPath)
  {
   //TODO:cant run without doing rectification yet
   if(NEED_RECTIFICATION)
   {

    double const maxZ = 1.0e4;
     
    ofstream ostr;
    string pointCloudFilename = m_params.resDir+filenameNoPath;
    ostr.open(pointCloudFilename.c_str(), ios::out);    

    for (int y = 0; y < m_points.rows; ++y) {
      for (int x = 0; x < m_points.cols; ++x) {
       
        Vec3f const& point = m_points.at<Vec3f>(y, x);
        //if (!(fabs(point[2] - maxZ) < FLT_EPSILON || fabs(point[2]) > maxZ)) {
        if (fabs(point[2]-maxZ)> FLT_EPSILON && point[2]<maxZ && point[2]>0){
          ostr << point[0] << " " << point[1] << " " << point[2] << endl;
        }
      }
    }
    ostr.close();
   }

  }

  //modifies the Q matrix with the rotation  and translation
  void CvStereoSGBMProcessor::changeCoordinates(cv::Mat const &rotationMat, cv::Mat const &vectorMat)
  {
   //TODO:cant do without doing rectification yet
   if(NEED_RECTIFICATION)
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
    gemm(extRotationMat, m_QOrig, 1.0,  m_Q, 0.0, m_Q, flags);
   
   }

  }

  void
  CvStereoSGBMProcessor::processImagePair(IplImage *refImg, IplImage *matchImg)
  {
    Mat refMatOrig = refImg;
    Mat matchMatOrig = matchImg;
    Mat tileDisparity;
    Rect imageRoi, tileRoi;
    Mat fullDisMat;
    Mat refMat, matchMat;

    // downsample
    resize(refMatOrig, refMat, smallSize, 0, 0, INTER_LINEAR );
    resize(matchMatOrig, matchMat, smallSize, 0, 0, INTER_LINEAR );

    if(NEED_RECTIFICATION)
    {
      //image rectification before stereo correlation
      remap(refMat, m_refImgRect, m_refRMap[0], m_refRMap[1], CV_INTER_LINEAR);
      remap(matchMat, m_matchImgRect, m_matchRMap[0], m_matchRMap[1], CV_INTER_LINEAR);
    }
    else
    {
    m_refImgRect = refMat.clone();
    m_matchImgRect = matchMat.clone();
    }    

    IplImage refImgRect;
    IplImage matchImgRect;

    refImgRect = m_refImgRect;
    matchImgRect = m_matchImgRect;

    fullDisMat = m_imageDisparity;

    //create tiles
    Tiling tiler;
    initializeTiling(tiler);
    tiler.setUpReferenceTiles(&refImgRect);

    //process tiles
    for(vector<CvRect>::iterator it = tiler.tiles.begin(); it != tiler.tiles.end(); ++it)
       {
       Mat refTile( m_refImgRect, *it );
       Mat matchTile( m_matchImgRect, *it );
       tileDisparity.create(it->height, it->width, CV_16S);

       //process tile
       sgbm(refTile,matchTile,tileDisparity);

       //create roi's for copying part of the tile without the overlap
       imageRoi = getTileWithoutOverlap(*it, tiler.tileParams.xOverlap, tiler.tileParams.yOverlap);
       tileRoi.width = imageRoi.width;
       tileRoi.height = imageRoi.height;
       tileRoi.x = imageRoi.x - it->x;
       tileRoi.y = imageRoi.y - it->y;

       //copy data using roi's
       Mat tileRoiMat(tileDisparity, tileRoi);
       Mat imageRoiMat(fullDisMat, imageRoi);
       tileRoiMat.copyTo(imageRoiMat);
       }

    // hack to get rid of small-scale noise
    IplImage disp = fullDisMat;
    
    cvErode(&disp, &disp, NULL, 2);
    cvDilate(&disp, &disp, NULL, 2);
    Mat disMat = &disp;

    //TODO:cant do this without doing rectification yet
    if(NEED_RECTIFICATION)
      reprojectImageTo3D(disMat, m_points, m_Q, 1);
  }

  //save the normalized disparity for display and debug  
  void
  CvStereoSGBMProcessor::saveDisparity(string const& filenameNoPath)
  {
    string disparityFilename = m_params.resDir+filenameNoPath;
    normalize(m_imageDisparity, m_imageDisparityNormalized, 0, 256, CV_MINMAX);
    imwrite( disparityFilename.c_str(), m_imageDisparityNormalized);
  }
  
cv::Mat CvStereoSGBMProcessor::GetRectifiedModelImage()
{
   return m_refImgRect;
}

cv::Mat CvStereoSGBMProcessor::GetRectifiedMatchImage()
{
   return m_matchImgRect;
}

// function to compute overlap and set up tile size
void CvStereoSGBMProcessor::initializeTiling(Tiling &tiling)
{
    // check for default no tiling
    if(m_params.tileWidth <= 0 || m_params.tileHeight <= 0)
      {
      m_params.tileWidth = smallSize.width;
      m_params.tileHeight = smallSize.height;
      }

    // set tile size from params
    tiling.tileParams.tileWidth = m_params.tileWidth;
    tiling.tileParams.tileHeight = m_params.tileHeight;

    // force tiles to be of even size
    if( tiling.tileParams.tileWidth%2 == 1 ) 
       tiling.tileParams.tileWidth++;
    if( tiling.tileParams.tileHeight%2 == 1 ) 
       tiling.tileParams.tileHeight++;

    // calculate overlap
    int minDisparity = m_params.minDisparity;
    int maxDisparity = m_params.minDisparity + m_params.numberOfDisparities;
    int largerDisparity = max( abs(minDisparity), abs(maxDisparity) );

    int xOverlap = largerDisparity * 2 + m_params.sadWindowSize + 1;
    int yOverlap = m_params.sadWindowSize + 1;

    // check tile size
    if(tiling.tileParams.tileWidth <= xOverlap)
      {
         cout<<"WARNING: tileWidth <= required xOverlap: "<<xOverlap<<endl<<"Setting tileWidth to 2*xOverlap"<<endl;
         tiling.tileParams.tileWidth = 2*xOverlap;
      }
    if(tiling.tileParams.tileHeight <= yOverlap)
      {
         cout<<"WARNING: tileHeight <= required yOverlap: "<<yOverlap<<endl<<"Setting tileWidth to 2*yOverlap"<<endl;
         tiling.tileParams.tileHeight = 2*yOverlap;
      }
    
    // set overlap
    tiling.tileParams.xOverlap = xOverlap;
    tiling.tileParams.yOverlap = yOverlap;

}



Rect CvStereoSGBMProcessor::getTileWithoutOverlap(Rect tile, int xOverlap, int yOverlap)
{
   Rect roi;

   if( tile.x == 0 ){
      roi.x = 0;
      roi.width = tile.width;
      }
   else{
      roi.x = tile.x + xOverlap/2;
      roi.width = tile.width - xOverlap/2;
      }
   if( tile.y == 0 ){
      roi.y = 0;
      roi.height = tile.height;
      }
   else{
      roi.y = tile.y + yOverlap/2;
      roi.height = tile.height - yOverlap/2;
      }

   return roi;
}

