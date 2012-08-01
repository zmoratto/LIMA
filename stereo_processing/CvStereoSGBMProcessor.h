// -*- C++ -*-

//#include "knStereo_Export.h"
//#include "StereoParameters.h"

// OpenCV includes
#include <opencv/cv.h>

#include <iosfwd>
#include <string>

// Tiling includes
#include "../common/Tiling.h"

struct CvStereoSGBMProcessorParameters{
  int preFilterCap;
  int sadWindowSize;
  int P1;
  int P2;
  int minDisparity;
  int numberOfDisparities;
  int uniquenessRatio;
  std::string modelImageFilename;
  std::string resDir;
  bool needRectification;
  float scaleFactor;
  int tileWidth;
  int tileHeight;
};

class CvStereoSGBMProcessor
{
 public:
  typedef cv::Mat PointImage;
  
  CvStereoSGBMProcessor(CvStereoSGBMProcessorParameters const& params);
  ~CvStereoSGBMProcessor();
  
  void processImagePair(IplImage *refImg, IplImage *matchImg);
  //void processImagePair(kn::Image const& refImg, kn::Image const& matchImg);
  
  cv::Mat getDisparityMap(){return m_imageDisparity;};
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
  void initializeTiling(Tiling &tiling);
  cv::Rect getTileWithoutOverlap(cv::Rect tile, int xOverlap, int yOverlap); 
 
  CvStereoSGBMProcessorParameters m_params;
  cv::StereoSGBM sgbm;
  
  cv::Mat m_refRMap[2];
  cv::Mat m_matchRMap[2];
  cv::Mat m_Q;
  cv::Mat m_QOrig;
  
  cv::Mat m_refImgRect;
  cv::Mat m_matchImgRect;
  
  cv::Mat m_imageDisparity;
  cv::Mat m_imageDisparityNormalized;
  
  cv::Size smallSize;

  PointImage m_points;

  bool NEED_RECTIFICATION;
};
