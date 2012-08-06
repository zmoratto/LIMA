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

#include "../SFM.h"
#include "../../stereo_processing/CvStereoBMProcessor.h"
#include "../../string_util.h"

using namespace std;

int ReadStereoConfigFile(string stereoConfigFilename, CvStereoBMProcessorParameters *thisStereoParams);

int main(int argc, char** argv)
{
	string modelList = "modelList.txt", matchList = "matchList.txt";
	char sfmConfigFile[] = "../sfm_config.txt";
	vector<string> modelImageFilenames;
	vector<string> matchImageFilenames;
	IplImage *modelImage = NULL, *matchImage = NULL;
	IplImage *image = NULL;
	IplImage im;
	cv::Mat rModelImage;
	int frameIndex, index;
	ifstream fin1, fin2;
	ofstream fout;
	string tempString;
	cv::Mat rotationMat=cv::Mat::eye(3,3, CV_64F);
	cv::Mat translationMat=cv::Mat::zeros(3,1, CV_64F);
	double const maxZ = 1.0e4;
		  
	//UBER Hack be careful - START
	float PI = 3.1415;
	float xAngleRad = 65*PI/180;
	//rotation.at<double>(0,0) = 1.0;
	rotationMat.at<double>(1,1) = cos(xAngleRad);
	rotationMat.at<double>(1,2) = -sin(xAngleRad);
	rotationMat.at<double>(2,1) = sin(xAngleRad);
	rotationMat.at<double>(2,2) = cos(xAngleRad);
	//UBER Hack be careful - END

	//Read Stereo Config File
	CvStereoBMProcessorParameters thisStereoParams;
	int configReadError = ReadStereoConfigFile(string("../stereo_settingsBM.txt"), &thisStereoParams);
	if( !configReadError )
	{
		cout << "Stereo Config File Read Error" << endl;
		return -1;
	}

	//Read Stereo Filenames
	fin1.open(modelList.c_str());
	fin2.open(matchList.c_str());
	while(fin1.good())
	{
		fin1 >> tempString;
		modelImageFilenames.push_back(tempString);
		fin2 >> tempString;
		matchImageFilenames.push_back(tempString);
	}
	fin1.close();
	fin2.close();

	//Set Up SFM
	SFM sfmTest(sfmConfigFile);
	image = cvLoadImage(modelImageFilenames[0].c_str(), 0);
	sfmTest.setUpSFM(NULL, image);

	//Set up Stereo
	CvStereoBMProcessor *thisStereo;
	thisStereo = new CvStereoBMProcessor(thisStereoParams);
	thisStereo->changeCoordinates(rotationMat, translationMat);


	for(frameIndex = sfmTest.configParams.firstFrame; frameIndex <= sfmTest.configParams.lastFrame; frameIndex++)
	{
		cout << endl << endl <<"Frame " << frameIndex << endl;
		//Load Stereo Images
		modelImage = cvLoadImage( modelImageFilenames[frameIndex].c_str(), 0 );
		matchImage = cvLoadImage( matchImageFilenames[frameIndex].c_str(), 0 );

		thisStereo->processImagePair(modelImage, matchImage);
		rModelImage = thisStereo->GetRectifiedModelImage();

		cv::Mat m_points = thisStereo->pointImage();

		//fout.open("results/pc.txt");
		//Save Point Clouds	 
		for (int y = 0; y < m_points.rows; ++y)
		{
			for (int x = 0; x < m_points.cols; ++x)
			{
				cv::Vec3f const& point = m_points.at<cv::Vec3f>(y, x);
				if (fabs(point[2]-maxZ)> FLT_EPSILON && point[2]<maxZ && point[2]>0)
				{
					index = y*sfmTest.imageWidth + x;
					sfmTest.pose.curr_x[index] = point[0];
					sfmTest.pose.curr_y[index] = point[1];
					sfmTest.pose.curr_z[index] = point[2];
					//fout << point[0] << " " << point[1] << " " << point[2] << endl;
				}
			}
		}
		//fout.close();

		im = rModelImage;

		sfmTest.process(&im, frameIndex);

		cvReleaseImage(&modelImage);
		cvReleaseImage(&matchImage);
		modelImage = NULL;
		matchImage = NULL;
	}
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
    sline9 >> identifier >> thisStereoParams->calibrationFilename;
    cout<<thisStereoParams->calibrationFilename<<endl;

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline10; 
    sline10<<line;
    sline10 >> identifier;
    if( sline10 >> val ){
       cout<<val<<endl;
       thisStereoParams->tileWidth=val;
       }
    else{
       cout<<"TILE_Width: "<<-1<<endl;
       thisStereoParams->tileWidth=-1;
       }

    getline (configFile,line);
    cout<<line<<endl;
    stringstream sline11; 
    sline11<<line;
    sline11 >> identifier;
    if( sline11 >> val ){
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

