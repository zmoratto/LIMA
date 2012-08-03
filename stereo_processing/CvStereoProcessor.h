// -*- C++ -*-

//#include "knStereo_Export.h"
//#include "StereoParameters.h"

// OpenCV includes
#include <opencv/cv.h>

#include <iosfwd>
#include <string>

struct CvStereoProcessorParameters{
  int preFilterSize;
  int preFilterCap;
  int sadWindowSize;
  int minDisparity;
  int numberOfDisparities;
  int textureThreshold;
  int uniquenessRatio;
  std::string modelImageFilename;
  std::string resDir;
  std::string calibrationFilename;
};

class CvStereoProcessor
{
 public:
  typedef cv::Mat PointImage;
  
  CvStereoProcessor(CvStereoProcessorParameters const& params);
  ~CvStereoProcessor();
  
  float samplingFactor;
  void processImagePair(IplImage *refImg, IplImage *matchImg);
  //void processImagePair(kn::Image const& refImg, kn::Image const& matchImg);
  
  CvMat* getDisparityMap(){return m_imageDisparity;};
  void saveDisparity(std::string const& filenameNoPath);
  
  void changeCoordinates(cv::Mat const& rotationMat, cv::Mat const& vectorMat);
  
  PointImage const& pointImage() const {
    return m_points;
  }
  

  void savePoints(std::string const& filenameNoPath);
  cv::Mat GetRectifiedModelImage();
  cv::Mat GetRectifiedMatchImage();
  
 protected:
  //static void iplImageFromKnImage(IplImage& ipl, Image const& img);
  
  CvStereoProcessorParameters m_params;
  CvStereoBMState * m_bMState;
  
  cv::Mat m_refRMap[2];
  cv::Mat m_matchRMap[2];
  cv::Mat m_Q;
  cv::Mat m_QOrig;
  
  cv::Mat m_refImgRect;
  cv::Mat m_matchImgRect;
  
  CvMat * m_imageDisparity;
  CvMat * m_imageDisparityNormalized;
  
  PointImage m_points;
};
