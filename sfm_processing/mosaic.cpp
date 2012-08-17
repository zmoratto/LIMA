// mosaic.cpp : mosaicing images and point cloud using tiles

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>

// OpenCV includes
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <opencv/cvaux.h>

// Program includes
#include "../camera_models/CAHV-to-pinhole.h"
#include "../common/opencv_pds.h"
#include "mosaic.h"

using namespace std;


// Erases a file extension if one exists and returns the base string
string
prefix_from_filename(string const& filename)
{
  string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
      result.erase(index, result.size());
  return result;
}

// Finds the file extension
string
extension_from_filename(string const& filename)
{
  string result = filename;
  int index = result.rfind(".");
  if (index != -1) 
    result.erase(0, index);
  return result;
}

// function to read the settings
struct mosaicSettings ReadMosaicSettings(const string &settingsFilename) 
{
  mosaicSettings settings;
  ifstream file;
  string line;
 
  //cout << "ReadCamFile(): camFilename = " << filename << endl;
  file.open(settingsFilename.c_str());
  if (file!=NULL){

    std::getline(file, line);
    int lineIndex = 0;
    int rows, cols;
    while (!file.eof()){
      if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines
      
      stringstream lineStream(line);
      string tmp, tmp1;
      
      if (lineIndex == 0){
	lineStream >> tmp >> settings.pixelSize;
        //cout<<tmp<<tmp1<<" "<<settings.pixelSize<<endl;
      }
     
      if (lineIndex == 1){
	lineStream >> tmp >> settings.tileWidth;
        //cout<<tmp<<", "<<settings.tileWidth<<endl;
      }
      if (lineIndex == 2){
	lineStream >> tmp >> settings.tileHeight;
        //cout<<tmp<<", "<<settings.tileHeight<<endl;
      }
      if (lineIndex == 3){
	lineStream >> tmp >> settings.scaleFactor;
        //cout<<tmp<<", "<<settings.scaleFactor<<endl;
      }
      if (lineIndex == 4){
	lineStream >> tmp >> settings.makeTileMosaic;
        //cout<<tmp<<", "<<settings.makeTileMosaic<<endl;
      }
      if (lineIndex == 5){
	lineStream >> tmp >> settings.makeImageMosaic;
        //cout<<tmp<<", "<<settings.makeImageMosaic<<endl;
      }
      lineIndex++;
      
      }
      std::getline(file, line);
    }
    file.close(); 
  }
  else{
    cout<<"WARNING! could not open settings file. use default values"<<endl;
    settings.pixelSize = 0.000012;
    settings.tileWidth = 512;
    settings.tileHeight = 512;
    settings.scaleFactor = 0.25;
    settings.makeTileMosaic = 1;
    settings.makeImageMosaic = 1;
  }
  return settings;
}

void PrintSettings(struct mosaicSettings settings)
{
  cout<<"pixelSize="<<settings.pixelSize<<endl;
  cout<<"tileWidth="<<settings.tileWidth<<endl;
  cout<<"tileHeight="<<settings.tileHeight<<endl;
  cout<<"scaleFactor="<<settings.scaleFactor<<endl;
  cout<<"makeTileMosaic="<<settings.makeTileMosaic<<endl;
  cout<<"makeImageMosaic="<<settings.makeImageMosaic<<endl;
}

mosaicProcessor::mosaicProcessor(const mosaicSettings &params)
{
   pixelSize = params.pixelSize; 
   tileWidth = params.tileWidth;
   tileHeight = params.tileHeight;
   scaleFactor = params.scaleFactor;
}

// function to determine which images overlap the tile
std::vector<int> mosaicProcessor::makeOverlapListFromBB(struct BBox tileBox, std::vector<BBox> imageBox)
{
  cout<<"NUM_BBoxes="<<imageBox.size()<<endl;

  std::vector<int> overlapIndices;
  overlapIndices.resize(imageBox.size());
  for (int i = 0; i< imageBox.size(); i++){
 
    if ((tileBox.maxPan < imageBox[i].minPan)||(imageBox[i].maxPan < tileBox.minPan)
	||(tileBox.maxTilt < imageBox[i].minTilt)||(imageBox[i].maxTilt < tileBox.minTilt)){
       overlapIndices[i] = 0;
    }
    else{
       overlapIndices[i] = 1;
       cout<<"****i:"<<i<<" minPan, maxPan, minTilt, maxTile"<<endl;
       cout<<"tile :"<<tileBox.minPan<<", "<<tileBox.maxPan<<", "<<tileBox.minTilt<<", "<<tileBox.maxTilt<<endl;
       cout<<"image:"<<imageBox[i].minPan<<", "<<imageBox[i].maxPan<<", "<<imageBox[i].minTilt<<", "<<imageBox[i].maxTilt<<endl;
    }
  }
  return overlapIndices;

}

// function to display the tile mosaic for debug
void mosaicProcessor::displayTileMosaic(string resultsDir)
{
  //check the output mosaick
  int firstHorTile = 0;
  int lastHorTile = 40;
  int firstVerTile = 0;
  int lastVerTile = 8;
  int newTileWidth = (int)floor(tileWidth*scaleFactor);//128; 
  int newTileHeight = (int) floor(tileHeight*scaleFactor);//128;
 
  int mosaicWidth = newTileWidth*(lastHorTile-firstHorTile+1);
  int mosaicHeight = newTileHeight*(lastVerTile-firstVerTile+1);
  cv::Mat mosaicMat(cv::Size(mosaicWidth, mosaicHeight),CV_8UC1);
  for (int i = firstHorTile; i <= lastHorTile; i++){
    for (int j = firstVerTile; j <= lastVerTile; j++){
     	string tileFilename;
	ostringstream ss_1, ss_2;
	ss_1 << i;
	ss_2 << j;
	tileFilename = resultsDir+"tile_"+ss_1.str()+"_"+ss_2.str()+".jpg";
	cout<<tileFilename<<endl;
     
	cv::Mat tileMat = cv::imread(tileFilename, 0);
	if (tileMat.data==NULL){
	  cout<<"i = "<<i<<", j = "<<j<<endl;
        }
        else{
	  cv::Mat smallTileMat(cv::Size(newTileWidth, newTileHeight),CV_8UC1);
	  cv::resize(tileMat, smallTileMat, smallTileMat.size(), 0, 0);
	  
	  int colStart = (i-firstHorTile)*newTileWidth;
          int rowStart = (j-firstVerTile)*newTileHeight;

          cout<<"colStart="<<colStart<<", rowStart="<<rowStart<<", newTileWidth="<<newTileWidth<<", newTileHeight="<<newTileHeight<<", mosaicWidth = "<<mosaicWidth<<", mosaicHeight="<<mosaicHeight<<endl;
	  cv::Mat mosaic_roi = mosaicMat(cv::Rect(colStart, rowStart, newTileWidth, newTileHeight));
          smallTileMat.copyTo(mosaic_roi);
       }
    }
  }
  string mosaicFilename = resultsDir+"tile_mosaic.jpg";
  cv::imwrite(mosaicFilename, mosaicMat);
}

//reads a list file and creates a vector of image pairs
//vector<struct imgPair> ReadListFile(string listFilename);
//reads a list file and creates a vector of image pairs
vector<ImagePairNames> mosaicProcessor::ReadListFile(string listFilename)
{
  ifstream file;
  string line;
  vector<ImagePairNames> imgPairVector;

  char *leftImgFilename;
  char *rightImgFilename;
  cout << endl;
  cout << "ReadListFile(): listFilename = " << listFilename << endl;
  file.open(listFilename.c_str());

  // Paranoid call to make sure the size really is zero... -- LJE
  imgPairVector.clear();

  std::getline(file, line);
  while (!file.eof())
  {
    if (!line.empty() && line.at(0) != '#') // skip comment & empty lines
    {
      stringstream lineStream(line);
      ImagePairNames thisImgPair;

      lineStream >> thisImgPair.left >> thisImgPair.right;
      if (lineStream)
	imgPairVector.push_back(thisImgPair);
    }
    std::getline(file, line);
  }
  file.close(); 
  cout << "ReadListFile(): imgPairVector.size() = "
       << imgPairVector.size() << endl;

  return imgPairVector;
}



std::vector<point3D> mosaicProcessor::readXYZFile(const string &inputFilename, int rows, int cols)
{
  string line;
  float x, y, z;
  ifstream myfile (inputFilename.c_str());
  std::vector<point3D> pc;

  pc.resize(cols*rows+1);
  point3D thisPC;

  if (myfile.is_open()){
    int index = 0;
    while ( myfile.good() ){
      getline (myfile,line);
      stringstream sline; 
      sline<<line;
      sline>> x >> y>> z;
   
      pc[index].x = x;
      pc[index].y = y;
      pc[index].z = z;
      
      index++;
    }
    cout<<"numPoints ="<<pc.size()<<endl;
    myfile.close();
  }
  else cout << "Unable to open file "<<inputFilename<<endl; 
  return pc;
}

void mosaicProcessor::saveXYZFile(const string &inputFilename, std::vector<point3D> &points, unsigned char *uchar_tile)
{
  ofstream myfile;
  myfile.open (inputFilename.c_str());
  for (int i=0; i < points.size(); i++){
    myfile <<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<" "<<(int)uchar_tile[i]<<" "<<(int)uchar_tile[i]<<" "<<(int)uchar_tile[i]<<endl;
  }
  myfile.close();
}

// function to read the cam params from a .cahv file
void mosaicProcessor::ReadCamParamsFromFile(char *filename, double*c, double*a, double*h, double*v, int* r_rows, int* r_cols)
{

  ifstream file;
  string line;
 
  //cout << "ReadCamFile(): camFilename = " << filename << endl;
  file.open(filename);
  if(!file.good()){
    cout<<"Error: file not opened: "<<filename<<endl;
    return;
      }

  std::getline(file, line);
  int lineIndex = 0;
  int rows, cols;
  while (!file.eof()){
    if (!line.empty() && line.at(0) != '#'){ // skip comment & empty lines
      stringstream lineStream(line);
      string tmp, tmp1;
      
      if (lineIndex == 0){
	lineStream >> tmp >> tmp1 >> c[0] >> c[1] >> c[2];
        //cout<<tmp<<tmp1<<c[0]<<", "<<c[1]<<", "<<c[2]<<endl;
      }
      if (lineIndex == 1){
	lineStream >> tmp >> tmp1 >> a[0] >> a[1] >> a[2];
        //cout<<tmp<<tmp1<<a[0]<<", "<<a[1]<<", "<<a[2]<<endl;
      }
      if (lineIndex == 2){
	lineStream >> tmp >> tmp1 >> h[0] >> h[1] >> h[2];
      }
      if (lineIndex == 3){
	lineStream >> tmp >> tmp1 >> v[0] >> v[1] >> v[2];
      }
      if (lineIndex == 4){
	lineStream >> tmp >> tmp1 >> cols;
      }
      if (lineIndex == 5){
	lineStream >> tmp >> tmp1 >> rows;
      }
      lineIndex++; 
    }
    std::getline(file, line);
  }
  file.close(); 
  *r_rows = rows;
  *r_cols = cols;
 
}

// function to take a quaternion and compute the 3x3 rotation matrix
void mosaicProcessor::quaternionToRotation( double quaternion[4], double matrix[3][3] )
{
//normalize
double mag = sqrt(quaternion[0]*quaternion[0]+quaternion[1]*quaternion[1]+quaternion[2]*quaternion[2]+quaternion[3]*quaternion[3]);

double w = quaternion[0]/mag, x = quaternion[1]/mag, y = quaternion[2]/mag, z = quaternion[3]/mag;

double xx = x*x, xy = x*y, xz = x*z, xw = x*w;
double yy = y*y, yz = y*z, yw = y*w;
double zz = z*z, zw = z*w;

matrix[0][0] = 1 - 2 * ( yy + zz );
matrix[0][1] = 2 * ( xy - zw );
matrix[0][2] = 2 * ( xz + yw );

matrix[1][0] = 2 * ( xy + zw );
matrix[1][1] = 1 - 2 * ( xx + zz );
matrix[1][2] = 2 * ( yz - xw );

matrix[2][0] = 2 * ( xz - yw );
matrix[2][1] = 2 * ( yz + xw );
matrix[2][2] = 1 - 2 * ( xx + yy );

}

//rotates the input vector by the quaternion
void mosaicProcessor::rotateWithQuaternion( double input[3], double output[3], double quaternion[4] )
{

double M[3][3];
quaternionToRotation( quaternion, M );

double a0 = M[0][0]*input[0] + M[0][1]*input[1] + M[0][2]*input[2];
double a1 = M[1][0]*input[0] + M[1][1]*input[1] + M[1][2]*input[2];
double a2 = M[2][0]*input[0] + M[2][1]*input[1] + M[2][2]*input[2];

output[0] = a0;
output[1] = a1;
output[2] = a2;

}

// function to compute the dot product of 3d vectors
inline double mosaicProcessor::dotprod3( double A[3], double B[3] ){
return (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
}

// function to get the x and y positions in an image from the 3d point p
void mosaicProcessor::getXYFromCAHV( double c[3], double a[3], double h[3], double v[3], double &x, double &y, double p[3] )
{
   double diff[3];
   for(int k=0;k<3;++k)
      diff[k]=p[k]-c[k];

x = dotprod3(diff,h)/dotprod3(diff,a);
y = dotprod3(diff,v)/dotprod3(diff,a);

}

// function to read cam params from a PDS file
void mosaicProcessor::ReadCamParamsFromPDSFile(char *filename, double*c, double*a, double*h, double*v, int* r_rows, int* r_cols)
{

  // read cahv
  float f_c[3],f_a[3],f_h[3],f_v[3],f_o[3],f_r[3];
  bool gotOR;
  if (scan_file_for_cahv_params(filename, f_c, f_a, f_h, f_v, f_o, f_r, gotOR)) {
    cout<<"Could not open "<<filename<<endl;
    return;
    } 

  for (int k = 0; k < 3; k++){
    c[k]=(double)f_c[k];
    a[k]=(double)f_a[k];
    h[k]=(double)f_h[k];
    v[k]=(double)f_v[k];
    }

  long l_rows, l_cols;
  scan_file_for_image_size( filename, l_rows, l_cols );
  *r_rows = (int) l_rows;
  *r_cols = (int) l_cols;
 
}


// function to display the image mosaic corrected using the quaternion
void mosaicProcessor::displayCorrectedImageMosaic(string resultsDir, string dataDir, vector<ImagePairNames> imgPairVector, 
                                 vector<BBox> BBoxArray, struct BBox mosaicBBox,  float radPerPixX, float radPerPixY)
{

    float minRadTilt = mosaicBBox.minTilt;
    float minRadPan = mosaicBBox.minPan;
    float maxRadTilt = mosaicBBox.maxTilt;
    float maxRadPan = mosaicBBox.maxPan;
   
    double c[3],a[3],h[3],v[3];
    //double c_avg[3] = {0.435831, 0.0372165, -1.22579};
    int cols, rows;

    cout<<"radPerPixX="<<radPerPixX<<", radPerPixY="<<radPerPixY<<endl;

 
    float deltaRadTilt = maxRadTilt-minRadTilt;
    float deltaRadPan  = maxRadPan-minRadPan;
    cout<<"deltaRadTilt="<<deltaRadTilt<<", deltaRadPan="<<deltaRadPan<<endl;

    cout<<"maxRadTilt: "<<maxRadTilt<<" minRadTilt: "<<minRadTilt<<" radPerPixY: "<<radPerPixY<<endl;
    cout<<"maxRadPan: "<<maxRadPan<<" minRadPan: "<<minRadPan<<" radPerPixX: "<<radPerPixX<<endl;
    cout<<"scaleFactor: "<<scaleFactor<<endl;

    int mosaicHeight = (maxRadTilt-minRadTilt)/radPerPixY;
    int mosaicWidth = (maxRadPan-minRadPan)/radPerPixX;
    mosaicHeight = mosaicHeight*scaleFactor;
    mosaicWidth = mosaicWidth*scaleFactor;   

     cout<<"mosaicWidth="<<mosaicWidth<<", mosaicHeight="<<mosaicHeight<<endl;
     cv::Mat mosaicMat = cv::Mat::zeros(cv::Size(mosaicWidth, mosaicHeight),CV_8UC1);

    int numPairs = imgPairVector.size();
    cv::Mat tileMat;    

    // loop over pairs
    for (int k=0; k<numPairs; ++k)
         {
  double q[4];
         //get rotated camera parameters and read image
         string filenameExtension = extension_from_filename(imgPairVector[k].left);
         if ((filenameExtension.compare(string(".img")) == 0) || (filenameExtension.compare(string(".IMG")) == 0) ){
           ReadCamParamsFromPDSFile((char*)(dataDir + imgPairVector[k].left).c_str(), c, a, h, v, &cols, &rows);
      
           // get origin rotation 
           scan_file_for_origin_rotation((char*)(dataDir + imgPairVector[k].left).c_str(), q); 
           // rotate cahv parameters with quaternion
           rotateWithQuaternion(c,c,q);
           rotateWithQuaternion(a,a,q);
           rotateWithQuaternion(h,h,q);
           rotateWithQuaternion(v,v,q);
            
           // read image
           readPDSFileToMat(dataDir + imgPairVector[k].left, tileMat);

         }
         else{
           string camFilename = dataDir + prefix_from_filename(imgPairVector[k].left) + string(".cahv");
           ReadCamParamsFromFile((char*)camFilename.c_str(), c, a, h, v, &cols, &rows);
           string imgFilename = dataDir + imgPairVector[k].left;
           tileMat = cv::imread(imgFilename, 0);
         }

         if (tileMat.data==NULL){
            cout<<"File "<< dataDir + imgPairVector[k].left <<" not found"<<endl;  
            continue;
            }

         // determine rows and columns corresponding to BBox
         int col_start = (BBoxArray[k].minPan - minRadPan)/radPerPixX * scaleFactor;
         int col_end = (BBoxArray[k].maxPan - minRadPan)/radPerPixX * scaleFactor;

         int row_start = (BBoxArray[k].minTilt - minRadTilt)/radPerPixY * scaleFactor;
         int row_end = (BBoxArray[k].maxTilt - minRadTilt)/radPerPixY * scaleFactor;

         // iterate over the window in the image mosaic
         for (int row = row_start; row < row_end; ++row)
         for (int col = col_start; col < col_end; ++col)
            {

            // get i,j position
            int i=row % mosaicHeight;
            int j=col % mosaicWidth;
            if(j<0)
               j+=mosaicWidth;
            if(i<0)
               i+=mosaicHeight;

            // calculate the 3d position vector (p) from i,j
            double p[3];
            float radX = j*radPerPixX/scaleFactor+minRadPan;
            float radY = i*radPerPixY/scaleFactor+minRadTilt;
            p[0] = sin(-radX)*cos(radY);
            p[1] = cos(radX)*cos(radY);
            p[2] = sin(radY);

            double x, y;

            //check to make sure dont project 180 degrees off
            // dot product is less than 0 if angle is >90
            // this isn't needed if we use the boundary boxes
            if(dotprod3(p,a)<0)  
               continue;


            // get the position where the point is in the image
            double origin[3] = {0.0,0.0,0.0};
            getXYFromCAHV(origin, a, h, v, x, y, p);

            // if not in the image skip it
            if(x < 0 || x > cols || y < 0 || y > rows)
               continue;

            // linear interpolation to get value since x, y are floats
            float dx = x - floor(x);
            float dy = y - floor(y);
            uchar val = tileMat.at<uchar>(floor(y),floor(x))*(1-dy)*(1-dx) 
                       + tileMat.at<uchar>(floor(y),ceil(x))*(1-dy)*(dx)
                       + tileMat.at<uchar>(ceil(y),floor(x))*(dy)*(1-dx)
                       + tileMat.at<uchar>(ceil(y),ceil(x))*(dy)*(dx);

            // only copy value if the pixel wasn't in a black strip
            if(tileMat.at<uchar>(floor(y),floor(x))!=0 && tileMat.at<uchar>(floor(y),ceil(x))!=0
               && tileMat.at<uchar>(ceil(y),floor(x))!=0 && tileMat.at<uchar>(ceil(y),ceil(x))!=0 )
                  mosaicMat.at<uchar>(i,j) = val;

            }

         }

    string mosaicFilename = resultsDir+"img_mosaic.jpg";
    cv::imwrite(mosaicFilename, mosaicMat);
}


// original function to display the image mosaic
void mosaicProcessor::displayImageMosaic(string resultsDir, string dataDir, vector<ImagePairNames> imgPairVector, 
                        struct BBox mosaicBBox,  float radPerPixX, float radPerPixY)
{

    float minRadTilt = mosaicBBox.minTilt;
    float minRadPan = mosaicBBox.minPan;
    float maxRadTilt = mosaicBBox.maxTilt;
    float maxRadPan = mosaicBBox.maxPan;
   
    double c[3],a[3],h[3],v[3];
    int cols, rows;

    cout<<"radPerPixX="<<radPerPixX<<", radPerPixY="<<radPerPixY<<endl;

 
    float deltaRadTilt = maxRadTilt-minRadTilt;
    float deltaRadPan  = maxRadPan-minRadPan;
    cout<<"deltaRadTilt="<<deltaRadTilt<<", deltaRadPan="<<deltaRadPan<<endl;

    cout<<"maxRadTilt: "<<maxRadTilt<<" minRadTilt: "<<minRadTilt<<" radPerPixY: "<<radPerPixY<<endl;
    cout<<"maxRadPan: "<<maxRadPan<<" minRadPan: "<<minRadPan<<" radPerPixX: "<<radPerPixX<<endl;
    cout<<"scaleFactor: "<<scaleFactor<<endl;

    int mosaicHeight = (maxRadTilt-minRadTilt)/radPerPixY;
    int mosaicWidth = (maxRadPan-minRadPan)/radPerPixX;
    mosaicHeight = mosaicHeight*scaleFactor;
    mosaicWidth = mosaicWidth*scaleFactor;   

     cout<<"mosaicWidth="<<mosaicWidth<<", mosaicHeight="<<mosaicHeight<<endl;
     cv::Mat mosaicMat(cv::Size(mosaicWidth, mosaicHeight),CV_8UC1);

    int numPairs = imgPairVector.size();

    for (int i = 0; i < numPairs; i++){
       
      string camFilename = dataDir + prefix_from_filename(imgPairVector[i].left) + string(".cahv");
      cout<<camFilename<<endl;
      
      ReadCamParamsFromFile((char*)camFilename.c_str(), c, a, h, v, &cols, &rows);
      
      float rad_pan, rad_tilt;
      
      a2pan_tilt(a, rad_pan, rad_tilt);

      rad_pan = -rad_pan;
      rad_tilt = -rad_tilt;

      double optical_center[2];    
      double fov_X_rad;
      double fov_Y_rad;
      
      GetPinholeModel(c, a, h, v, cols, rows, pixelSize, optical_center, &fov_X_rad, &fov_Y_rad);
   
      string imgFilename = dataDir + imgPairVector[i].left;
      cv::Mat tileMat = cv::imread(imgFilename, 0);

      if ((tileMat.data==NULL)){
	cout<<"File "<<imgFilename<<" not found"<<endl;  
      }
      else{

	int newTileWidth = cols*scaleFactor;
        int newTileHeight = rows*scaleFactor;
        double new_optical_center[2];
        new_optical_center[0] = optical_center[0]*scaleFactor;
	new_optical_center[1] = optical_center[1]*scaleFactor;
        float new_rad_per_pix_X = radPerPixX/scaleFactor;
	float new_rad_per_pix_Y = radPerPixY/scaleFactor;

	cv::Mat smallTileMat(cv::Size(newTileWidth, newTileHeight),CV_8UC1);
	cv::resize(tileMat, smallTileMat, smallTileMat.size(), 0, 0);
      

	int colStart = (rad_pan-minRadPan)/(new_rad_per_pix_X) - new_optical_center[0];
     
	int rowStart = (rad_tilt-minRadTilt)/(new_rad_per_pix_Y) -  new_optical_center[1];
 
        if (rowStart<0) rowStart = 0;
        if (colStart<0) colStart = 0;

      
	cout<<"colStart="<<colStart<<", rowStart="<<rowStart<<", width="<<newTileWidth<<", height="<<newTileHeight<<", mosaicWidth="<<mosaicWidth<<", mosaicHeight="<<mosaicHeight<<endl;
	cout<<"rad_pan="<<rad_pan<<", rad_tilt="<<rad_tilt<<", minRadPan="<<minRadPan<<", maxRadPan="<<maxRadPan<<", minRadTilt="<<minRadTilt<<", maxRadTilt="<<maxRadTilt<<endl;
	if ((colStart+newTileWidth<mosaicWidth) && (rowStart+newTileHeight < mosaicHeight)){  
	    cv::Mat mosaic_roi = mosaicMat(cv::Rect(colStart, rowStart, newTileWidth, newTileHeight));
	    smallTileMat.copyTo(mosaic_roi);        
	}
	
      }
     
    }
    
    string mosaicFilename = resultsDir+"img_mosaic.jpg";
    cv::imwrite(mosaicFilename, mosaicMat);
}

// function that makes the tiles
void mosaicProcessor::makeTiles(string dataDir, string resultsDir, std::vector<ImagePairNames> imgPairVector, std::vector<BBox> BBoxArray, 
               struct BBox mosaicBBox, double radPerPixX, double radPerPixY)
{

  float minPanRad = (mosaicBBox.minPan);
  float maxPanRad = (mosaicBBox.maxPan);
  float minTiltRad = (mosaicBBox.minTilt);
  float maxTiltRad = (mosaicBBox.maxTilt);

  // check for default case in width
  if(tileWidth <= 0 || tileWidth*radPerPixX > maxPanRad - minPanRad )
      tileWidth = (maxPanRad - minPanRad)/radPerPixX;

  // check for default case in height
  if(tileHeight <= 0 || tileHeight*radPerPixY > maxTiltRad - minTiltRad )
      tileHeight = (maxTiltRad - minTiltRad)/radPerPixY;

  float tile_fov_x = tileWidth*radPerPixX;
  float tile_fov_y = tileHeight*radPerPixY;
  cout<<"tile_fov_x="<<tile_fov_x<<", tile_fov_y="<<tile_fov_y<<endl;

  int numHorTiles = ceil((maxPanRad-minPanRad)/tile_fov_x);
  int numVerTiles = ceil((maxTiltRad-minTiltRad)/tile_fov_y);
  cout<<"numHorTiles="<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;

  // loop over the tiles
  for (int i = 0; i < numHorTiles; i++){
    for (int j = 0; j < numVerTiles; j++){
    
      cout<<"i="<<i<<", j="<<j<<endl;

     // get the tile boundary box
     struct BBox tileBox;
     tileBox.minPan = (i*tileWidth)*radPerPixX + minPanRad;
     tileBox.maxPan = (i+1)*tileWidth*radPerPixX + minPanRad;
     tileBox.minTilt = j*tileHeight*radPerPixY + minTiltRad;
     tileBox.maxTilt = (j+1)*tileHeight*radPerPixY + minTiltRad;

     // get overlap list
     std::vector<int> overlapIndices = makeOverlapListFromBB(tileBox, BBoxArray);
  
     float *imgTile    = new float[tileWidth*tileHeight];
     unsigned char  *uchar_tile    = new unsigned char[tileWidth*tileHeight];
     int   *counter = new int[tileWidth*tileHeight]; 
     int   *xyzCounter = new int[tileWidth*tileHeight];
     std::vector<point3D> xyzTile;
     xyzTile.resize(tileWidth*tileHeight);
     
     //initialize the tile - START
     for (int l = 0; l < tileHeight; l++){//ver    
       for (int k = 0; k < tileWidth; k++){//hor 
	 int tileIndex = l*tileWidth+k;
	 imgTile[tileIndex] = 0.0;
	 uchar_tile[tileIndex] = 0;
	 counter[tileIndex] = 0;
	 xyzCounter[tileIndex] = 0;
	 xyzTile[tileIndex].x = 0.0;
	 xyzTile[tileIndex].y = 0.0;
	 xyzTile[tileIndex].z = 0.0;
       }
     }
     //initialize the tile - END

   
      //update tile - START
      int numPairs = imgPairVector.size();
      for (int p = 0; p < numPairs; p++){

	unsigned char *buffer;
	double fov_X, fov_Y; 
	int rows;
	int cols;
	double optical_center[2];

	if (overlapIndices[p] == 1){
	  
	  cout<<"p="<<p<<", "<<imgPairVector[p].left<<endl;

	  //load individual images - START
     IplImage *leftImage;
	  string imgFilename = dataDir + imgPairVector[p].left;
     string filenameExtension = extension_from_filename(imgPairVector[i].left);
     if ((filenameExtension.compare(string(".img")) == 0) || (filenameExtension.compare(string(".IMG")) == 0) ){
       readPDSFileToIplImage(imgFilename, leftImage);
     }
     else{
      leftImage = cvLoadImage((imgFilename).c_str());
     }

	  int t_rows = leftImage->height;
	  int t_cols = leftImage->width;
       
	  buffer = new unsigned char[t_rows*t_cols];
	  for (int y = 0; y < t_rows; y++){
	    for (int x = 0; x < t_cols; x++){
	      buffer[y*t_cols+x] = (unsigned char)(cvGet2D(leftImage, y, x).val[0]);
	    }
	  }
	  //buffer = (unsigned char*)(leftImage->imageData);
	  //load individual images - END
	  
	  //load individual point cloud files - START
	  string xyzFilename = resultsDir + string("point_") + prefix_from_filename(imgPairVector[p].left)+string(".txt"); 
	  cout<<"reading xyz file...";
          std::vector<point3D> points = readXYZFile(xyzFilename, t_rows, t_cols);
	  cout<<"done."<<endl;  
	  //load individual point cloud files - END

     //read cam params
	  double c[3],a[3],h[3],v[3];
     if ((filenameExtension.compare(string(".img")) == 0) || (filenameExtension.compare(string(".IMG")) == 0) ){
       ReadCamParamsFromPDSFile((char*)(dataDir + imgPairVector[p].left).c_str(), c, a, h, v, &cols, &rows);
       // get origin rotation quaternion
       double q[4];
       scan_file_for_origin_rotation((char*)(dataDir + imgPairVector[p].left).c_str(), q); 
       // rotate cahv parameters with quaternion
       rotateWithQuaternion(c,c,q);
       rotateWithQuaternion(a,a,q);
       rotateWithQuaternion(h,h,q);
       rotateWithQuaternion(v,v,q);
     }
     else{
       string camFilename = dataDir + prefix_from_filename(imgPairVector[p].left) + string(".cahv");
       ReadCamParamsFromFile((char*)camFilename.c_str(), c, a, h, v, &cols, &rows);
     }
	  
	  float cam_pan_rad, cam_tilt_rad;
	  a2pan_tilt(a, cam_pan_rad, cam_tilt_rad);
	  cam_pan_rad = -cam_pan_rad;
	  cam_tilt_rad = -cam_tilt_rad;
	  
	 
	  double fov_X_rad, fov_Y_rad; 
	  GetPinholeModel(c, a, h, v, cols, rows, pixelSize, optical_center, &fov_X_rad, &fov_Y_rad);

	  //fill in the tile	
	  for (int l = 0; l < tileHeight; l++){
	    for (int k = 0; k < tileWidth; k++){
	      
	      float curr_pan_rad  = ((i*tileWidth  + k)*radPerPixX + minPanRad);//absolute coordinates of the curr_pan
	      float curr_tilt_rad = ((j*tileHeight + l)*radPerPixY  + minTiltRad);//absolute coordinates of curr_tilt
	            
	      int pixVal;
	      struct point3D xyz;

	      int x, y;
	      //x, y are the row and column in the current image corresponding to curr_pan and curr_tilt
	      //x = (int)floor(optical_center[0] + (curr_pan_rad  - cam_pan_rad)/radPerPixX);
	      //y = (int)floor(optical_center[1] + (curr_tilt_rad - cam_tilt_rad)/radPerPixY);

         // calculate the 3d position vector (p) from current pan and tilt
         double p[3];
         p[0] = sin(-curr_pan_rad)*cos(curr_tilt_rad);
         p[1] = cos(curr_pan_rad)*cos(curr_tilt_rad);
         p[2] = sin(curr_tilt_rad);
         double x_d,y_d;

         // get the position where the point is in the image
         double origin[3] = {0.0,0.0,0.0};
         getXYFromCAHV(origin, a, h, v, x_d, y_d, p);
         x = (int)x_d;
         y = (int)y_d;

         // get imgIndex and tileIndex
         int imgIndex = y*cols+x;
         int tileIndex = l*tileWidth+k;
     
         //if the point lies inside the image 
         if ((x >=0) && (y >=0) && (x < cols) && (y <rows)){
         int pixVal = (int)(buffer[imgIndex]);

         //compute the distance the point is from the origin
         double dist = sqrt(points[imgIndex].x*points[imgIndex].x+points[imgIndex].y*points[imgIndex].y+points[imgIndex].z*points[imgIndex].z);

         // if the pixel is not 0
         if (pixVal!=0){
            imgTile[tileIndex]   = imgTile[tileIndex] + pixVal;
            counter[tileIndex]   = counter[tileIndex] + 1;

            // check to make sure the point wasnt (0, 0, 0)
            if(dist > FLT_EPSILON){
               xyzTile[tileIndex].x = xyzTile[tileIndex].x + p[0]*dist;
               xyzTile[tileIndex].y = xyzTile[tileIndex].y + p[1]*dist;
               xyzTile[tileIndex].z = xyzTile[tileIndex].z + p[2]*dist;
               xyzCounter[tileIndex] = xyzCounter[tileIndex] + 1;
               }
            }

         }
	    
       }
     } 
	  
	  cvReleaseImage(&leftImage);
	  delete[] buffer;
	  points.clear();
	  cout<<"done."<<endl;
	 }//if (overlapIndices[p] == 1){
	  
      }
      //update tile - END
      

      //normalize tile - START
      int numValidPix = 0;
      for (int l = 0; l < tileHeight; l++){
	for (int k = 0; k < tileWidth; k++){
	  int index = l*tileWidth+k;
          if (counter[index]>0){
	    uchar_tile[index] = (unsigned char) (imgTile[index]/counter[index]);
	    numValidPix++;
          }
       if(xyzCounter[index]>0){
         xyzTile[index].x = (xyzTile[index].x/xyzCounter[index]);
         xyzTile[index].y = (xyzTile[index].y/xyzCounter[index]);
         xyzTile[index].z = (xyzTile[index].z/xyzCounter[index]);
         }

	}
      }
      //normalize tile - END
      cout<<"done normalizing the tile"<<endl;

      //save tile - START
      if (numValidPix > 0){
	string tileFilename;
	ostringstream ss_1, ss_2;
	ss_1 << i;
	ss_2 << j;
	tileFilename = resultsDir+"tile_"+ss_1.str()+"_"+ss_2.str()+".jpg";
        string xyzTileFilename = resultsDir+"tile_"+ss_1.str()+"_"+ss_2.str()+".txt";

	cout<<tileFilename<<endl;
	IplImage *tileImage = cvCreateImage(cvSize(tileWidth, tileHeight), IPL_DEPTH_8U , 1);

	for (int k = 0; k < tileWidth; k++){
	  for (int l = 0; l < tileHeight; l++){
              int index = l*tileWidth+k;
	      ((uchar *)(tileImage->imageData +l*tileImage->widthStep))[k]=uchar_tile[index];
	  }
        }
        cvSaveImage(tileFilename.c_str(), tileImage);
	cvReleaseImage(&tileImage);
        saveXYZFile(xyzTileFilename, xyzTile, uchar_tile);
      }
      cout<<"done saving the tile"<<endl;
      //save tile - END

      delete[] uchar_tile;
      delete[] imgTile;
      delete[] counter;
      delete[] xyzCounter;
      xyzTile.clear();
  
    }
  }
}

//function to calculate boundary boxes for each image pair
void mosaicProcessor::calcBoundaryBoxes( string dataDir, vector<ImagePairNames> &imgPairVector, vector<BBox> &BBoxArray, BBox &mosaicBBox, float &radPerPixX, float &radPerPixY )
{

  int numPairs = imgPairVector.size();  
  BBoxArray.clear();

  float minPan = M_PI;
  float maxPan = -M_PI;
  float minTilt = M_PI/2.0;
  float maxTilt = -M_PI/2.0;
  
  //parameters characteristic to each image - START
  // we assume all images have the same fov, cols, rows...   
  double c[3],a[3],h[3],v[3],o[3],r[3];
  bool gotOR = false;
  int cols;
  int rows;
  double opticalCenter[2];    
  double radFovX;
  double radFovY;
  //parameters characteristic to each image - END
  
  //compute camera parameters
  cout<<"NUM_PAIRS="<<numPairs<<endl;

  for (int i = 0; i < numPairs; i++){
    
    //read the camera models - START
    string filenameExtension = extension_from_filename(imgPairVector[i].left);
    if ((filenameExtension.compare(string(".img")) == 0) || (filenameExtension.compare(string(".IMG")) == 0) ){
      ReadCamParamsFromPDSFile((char*)(dataDir + imgPairVector[i].left).c_str(), c, a, h, v, &cols, &rows);
      // get origin rotation quaternion
      double q[4];
      scan_file_for_origin_rotation((char*)(dataDir + imgPairVector[i].left).c_str(), q); 
      // rotate cahv parameters with quaternion
      rotateWithQuaternion(c,c,q);
      rotateWithQuaternion(a,a,q);
      rotateWithQuaternion(h,h,q);
      rotateWithQuaternion(v,v,q);
    }
    else{
      ReadCamParamsFromFile((char*)(dataDir + prefix_from_filename(imgPairVector[i].left) + string(".cahv")).c_str(), c, a, h, v, &cols, &rows);
    }
    //read the camera models - END

    // calculate boundary box - START
    float panRad, tiltRad;
 
    a2pan_tilt(a, panRad, tiltRad);
    panRad = -panRad;
    tiltRad = -tiltRad;
    cout<<"panRad="<<panRad<<", tiltRad="<<tiltRad<<endl;

    GetPinholeModel(c, a, h, v, cols, rows, pixelSize, opticalCenter, &radFovX, &radFovY);
    
    cout<<"rows: "<<rows<<endl;
    radPerPixX = radFovX/cols; 
    radPerPixY = radFovY/rows; 

    BBox thisBBox;   
    //FIXME: *1.5 was added because the boundary was too small
    thisBBox.minPan  = (panRad - radFovX/2*1.5);
    thisBBox.maxPan  = (panRad + radFovX/2*1.5);
    thisBBox.minTilt = (tiltRad - radFovY/2*1.5);
    thisBBox.maxTilt = (tiltRad + radFovY/2*1.5);
    BBoxArray.push_back(thisBBox);

    if (thisBBox.minTilt<minTilt){
      minTilt = thisBBox.minTilt;
    }
    if (thisBBox.maxTilt>maxTilt){
      maxTilt = thisBBox.maxTilt;
    }
    if (thisBBox.minPan<minPan){
      minPan = thisBBox.minPan;
    }
    if (thisBBox.maxPan>maxPan){
      maxPan = thisBBox.maxPan;
    }
       
  }

  minPan  = minPan  - radFovX/2;
  maxPan  = maxPan  + radFovX/2;
  minTilt = minTilt - radFovY/2;
  maxTilt = maxTilt + radFovY/2;

  cout<<"minPan="<<minPan<<", maxPan="<<maxPan<<", minTilt="<<minTilt<<", maxTilt="<<maxTilt<<endl;
  if (maxPan>M_PI){maxPan = M_PI;}
  if (minPan<-M_PI){minPan = -M_PI;} 
  if (maxTilt>M_PI/2){maxTilt = M_PI/2;}
  if (minTilt<-M_PI/2){minTilt = -M_PI/2;}
  mosaicBBox.minTilt = minTilt;
  mosaicBBox.maxTilt = maxTilt;
  mosaicBBox.minPan  = minPan;
  mosaicBBox.maxPan  = maxPan;

  // calculate boundary box - END


}



