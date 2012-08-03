#include <iostream>
#include <time.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

#include "../PoseEstimation.h"
#include "../FeatureExtraction.h"
#include "../SFM.h"
#include "../../stereo_processing/CvStereoBMProcessor.h"
#include "../../string_util.h"

using namespace std;

int ReadStereoConfigFile(string stereoConfigFilename, CvStereoBMProcessorParameters *thisStereoParams);


int main(int argc, char** argv)
{
	string modelImageFilename;
	string matchImageFilename;
	int verbose;
	string configFilename;
	string errorFilename;
	string resDir;
	string auxDir;

	// load images
	IplImage *modelImage, *matchImage;

	modelImage = cvLoadImage( modelImageFilename.c_str(), 0 );
	matchImage = cvLoadImage( matchImageFilename.c_str(), 0 );

	if(modelImage==NULL)
	{
		cerr<<"Error: image failed to load"<<endl;
	return -1;
	}
	if(matchImage==NULL)
	{
		cerr<<"Error: image failed to load"<<endl;
		return -1;
	}

	//run openCV stereo
	CvStereoBMProcessorParameters thisStereoParams;
	int configReadError = ReadStereoConfigFile(string("stereo_settingsBM.txt"), &thisStereoParams);
	 
	if( !configReadError )
		return -1;

	thisStereoParams.modelImageFilename = modelImageFilename;
	thisStereoParams.resDir = resDir;

	cv::Size imageSize;
	imageSize.width = modelImage->width;
	imageSize.height = modelImage->height;
	  
	CvStereoBMProcessor *thisStereo;  
	thisStereo = new CvStereoBMProcessor(thisStereoParams);


	cv::Mat rotationMat=cv::Mat::eye(3,3, CV_64F);
	cv::Mat translationMat=cv::Mat::zeros(3,1, CV_64F);
	  
	//UBER Hack be careful - START
	float PI = 3.1415;
	float xAngleRad = 65*PI/180;
	//rotation.at<double>(0,0) = 1.0;
	rotationMat.at<double>(1,1) = cos(xAngleRad);
	rotationMat.at<double>(1,2) = -sin(xAngleRad);
	rotationMat.at<double>(2,1) = sin(xAngleRad);
	rotationMat.at<double>(2,2) = cos(xAngleRad);
	//UBER Hack be careful - END

	thisStereo->changeCoordinates(rotationMat, translationMat);
	 
	thisStereo->processImagePair(modelImage, matchImage);
	IplImage rModelImage = thisStereo->GetRectifiedModelImage();

	cv::Mat m_points = thisStereo->pointImage();

	//SAVE POINT CLOUDS
	double const maxZ = 1.0e4;
		 
	for (int y = 0; y < m_points.rows; ++y)
	{
		for (int x = 0; x < m_points.cols; ++x)
		{
			cv::Vec3f const& point = m_points.at<cv::Vec3f>(y, x);
			//if (fabs(point[2]-maxZ)> FLT_EPSILON && point[2]<maxZ && point[2]>0)
			//	ostr << point[0] << " " << point[1] << " " << point[2] << endl;
		}
	}

	//END SAVE POINT CLOUDS
}

int ReadStereoConfigFile(string stereoConfigFilename, CvStereoBMProcessorParameters *thisStereoParams)
{
  ifstream configFile (stereoConfigFilename.c_str());
  std::string line;
  double val; 
  std::string identifier;
  std::string param;

  if (configFile.is_open()){ 
    
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline; 
    sline<<line;
    sline >> identifier >> val;
    cout<<val<<endl;
    thisStereoParams->preFilterSize=val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline1; 
    sline1<<line;
    sline1 >> identifier >> val;
      cout<<val<<endl;
    thisStereoParams->preFilterCap=val;
  
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline2; 
    sline2<<line;
    sline2 >> identifier >> val;
      cout<<val<<endl;
    thisStereoParams->sadWindowSize=val;
 
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline3; 
    sline3<<line;
    sline3 >> identifier >> val;
      cout<<val<<endl;
    thisStereoParams->minDisparity= val;
  
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline4; 
    sline4<<line;
    sline4 >> identifier >> val;
      cout<<val<<endl;
    thisStereoParams->numberOfDisparities= val;
  

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline5; 
    sline5<<line;
    sline5 >> identifier >> val;
   cout<<val<<endl;
    thisStereoParams->textureThreshold=val;
 
    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline6; 
    sline6<<line;
    sline6 >> identifier >> val;
      cout<<val<<endl;
    thisStereoParams->uniquenessRatio=val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline7; 
    sline7<<line;
    sline7 >> identifier >> param;
      cout<<param<<endl;
    if(param[0]=='N' || param[0]=='n')
      thisStereoParams->needRectification=false;
    else if(param[0]=='Y' || param[0]=='y')
      thisStereoParams->needRectification=true;
    else{
      cout << "Error reading stereo settings file"<<endl;
      configFile.close();
      return 0;
    }

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline8; 
    sline8<<line;
    sline8 >> identifier >> val;
    cout<<val<<endl;
    thisStereoParams->scaleFactor=val;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline9; 
    sline9<<line;
    sline9 >> identifier;
    if( sline9 >> val ){
       cout<<val<<endl;
       thisStereoParams->tileWidth=val;
       }
    else{
       cout<<"TILE_HEIGHT: "<<-1<<endl;
       thisStereoParams->tileWidth=-1;
       }

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline10; 
    sline10<<line;
    sline10 >> identifier;
    if( sline10 >> val ){
       cout<<val<<endl;
       thisStereoParams->tileHeight=val;
       }
    else{
       cout<<"TILE_HEIGHT: "<<-1<<endl;
       thisStereoParams->tileHeight=-1;
       }

    configFile.close();
    return 1;
  }
  else{
    cout << "Unable to open settings file"<<endl; 
    return 0;
  } 
}
