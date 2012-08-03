#include "../CvStereoBMProcessor.h"
#include "../../string_util.h"

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
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


using std::cout;
using std::endl;

bool verbose = true;


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


//reads the calibrartion file
int ReadCalibrationFile(string calibrationFilename, cv::Mat &camera_matrix_left, cv::Mat &camera_matrix_right,
                                                    cv::Mat &dist_coeffs_left, cv::Mat &dist_coeffs_right, 
                                                    cv::Mat &R, cv::Mat &T)
{
  ifstream calibrationFile (calibrationFilename.c_str());
  std::string line;
  double a00, a01,a02, a10, a11, a12, a20, a21, a22;

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
    //R = (Mat_<double>(3,3)<<a00, a01, a02, a10, a11, a12, a20, a21, a22);
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
    //cout << line << endl;
    sline5 << line;
    sline5 >> identifier >> a00 >> a01 >> a02;
    cout<<identifier<<" "<<a00<<" "<<a01<<" "<<a02<<endl;
    //T = (Mat_<double>(3,1)<<a00, a01, a02);  
    T.create(3, 1, CV_64F);
    T.at<double>(0,0)=a00;
    T.at<double>(1,0)=a01;
    T.at<double>(2,0)=a02;
    //T.at<double>(3,0)=a10;

    calibrationFile.close();
    return 1;
  }
  else {
    cout << "Unable to open file"; 
    return 0;
  }
}

int main( int argc, char *argv[] )
{

 std::string modelImageFilename;
 std::string matchImageFilename;
 
 int verbose;
 std::string configFilename;
 std::string errorFilename;
 std::string resDir;
 std::string auxDir;

 po::options_description general_options("Options");
 general_options.add_options()
    ("model-filename,a", po::value<std::string>(&modelImageFilename))
    ("match-filename,b", po::value<std::string>(&matchImageFilename))
    ("results-directory,r", 
		po::value<std::string>(&resDir)->default_value("results"), 
		"results directory.") // Currently a no-op until we implement it below.
    ("settings-filename,s", 
		po::value<std::string>(&configFilename)->default_value("stereo_settings.txt"), 
		"settings filename.")
    ("error-filename,e", 
		po::value<std::string>(&errorFilename)->default_value("stereo_errors.txt"), 
		"error output filename.")
    ("verbose,v", 
		po::value<int>(&verbose)->default_value(1), 
		"Verbosity level, zero emits no messages.")
    ("help,h", "Display this help message");
 
  po::options_description hidden_options("");

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  //p.add("model-filename", -1);

  std::ostringstream usage;
  usage << "Description: main code for image to image co-registration" << std::endl << std::endl;
  usage << general_options << std::endl;
 
  po::variables_map vm;
  try
    {
      po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
      po::notify( vm );
    }
  catch ( po::error &e )
   {
     std::cout << "An error occured while parsing command line arguments.\n";
     std::cout << "\t" << e.what() << "\n\n";
     std::cout << usage.str();
     return 1;
   }
  
  if( vm.count("help") )
    {
      std::cerr << usage.str() << std::endl;
      return 1;
    }

  if (!vm.count("model-filename")){
      std::cerr << "Error: Must specify at least one model image file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
  }   
  if (!vm.count("match-filename")){
      std::cerr << "Error: Must specify at least one match image file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
  }


  // load images
  IplImage *modelImage, *matchImage;

  modelImage = cvLoadImage( modelImageFilename.c_str(), 0 );
  matchImage = cvLoadImage( matchImageFilename.c_str(), 0 );

  if(modelImage==NULL){
   cerr<<"Error: image failed to load"<<endl;
   return -1;
	}
  if(matchImage==NULL){
   cerr<<"Error: image failed to load"<<endl;
   return -1;
   }

  string makeResultsDirCmd = "mkdir " + resDir;
  resDir = resDir + "/" + GetFilenameNoExt(GetFilenameNoPath(modelImageFilename)) + "_" + GetFilenameNoExt(GetFilenameNoPath(matchImageFilename));
  auxDir = resDir+"/aux";
  cout<<"aux_dir="<<auxDir<<endl;
  //create the results and auxiliary directories
  string makeResDirCmd = "mkdir " + resDir;
  string makeAuxDirCmd = "mkdir " + auxDir;
  system(makeResultsDirCmd.c_str());
  system(makeResDirCmd.c_str()); 
  system(makeAuxDirCmd.c_str());

 
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

  //get start time for debug - START
  timeval start, end;
  long utime, seconds, useconds;
  gettimeofday(&start, NULL);
  //get start time for debug - END

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

  //get end time and print for debug - START
  gettimeofday(&end, NULL);
  seconds = end.tv_sec - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;
  utime = ((seconds) * 1000000 + useconds);
  printf("Elapsed time: %f seconds\n", float(utime)/1000000.0);
  //get end time and print for debug - END

  //saves disparity map for debug - START
  string disparityFilenameNoPath = string("/disp_") + GetFilenameNoExt(GetFilenameNoPath(modelImageFilename)) + string(".jpg");
  cout<<disparityFilenameNoPath<<endl;
  thisStereo->saveDisparity(disparityFilenameNoPath);
  //saves disparity map for debug - END

  //saves the point clouds for debug - START 
  string pointCloudFilenameNoPath = string("/point_") + GetFilenameNoExt(GetFilenameNoPath(modelImageFilename)) + string(".txt");
  cout<<pointCloudFilenameNoPath<<endl;
  thisStereo->savePoints(pointCloudFilenameNoPath);
  //saves the point clouds for debug - END*/

  //copies the input rectified images for debug - START
  string outModelImage = resDir + "/" + GetFilenameNoExt(GetFilenameNoPath(modelImageFilename)) + string(".jpg");
  string outMatchImage = resDir + "/" + GetFilenameNoExt(GetFilenameNoPath(matchImageFilename))  +string(".jpg");
  IplImage rModelImage = thisStereo->GetRectifiedModelImage();
  IplImage rMatchImage = thisStereo->GetRectifiedMatchImage();
  cvSaveImage( outModelImage.c_str(), &rModelImage);
  cvSaveImage( outMatchImage.c_str(), &rMatchImage);
  //copies the input images for debug - END

  delete thisStereo;

  cvReleaseImage(&modelImage);
  cvReleaseImage(&matchImage);
}
