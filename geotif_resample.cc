// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>
#include <getopt.h>

//boost
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>

#include <stdint.h>
#include "geotif_resample.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;
/*
struct TileParams
{
  int verOffset;
  int horOffset;
  int padImgWidth;
  int padImgHeight;
  int numPyrLevels;
};

struct ResampleParams
{
  int imageType; //0 (DEM) or 1 (DRG), temporary solution
  int numPyrLevels; //negative numbers
  float pyrResampleFactor;
  int imageROILeft;
  int imageROITop;
  int imageROIRight;
  int imageROIBottom;
  int tileSize;
  int tilePaddingTop;
  int tilePaddingLeft;
  int tilePaddingRight;
  int tilePaddingBottom;
};
*/


static bool
CPUIsBigEndian()
{
  union { short number; char bytes[2]; } short_endian_check;

  short_endian_check.number = 1;

  return (char(short_endian_check.number) == short_endian_check.bytes[1]);
}

static uint16_t
SwapWord(unsigned short src)
{
  uint16_t *shortPtr = &src;
  uint8_t temp = ((unsigned char *) shortPtr)[1];

  ((uint8_t *) shortPtr)[1] = ((uint8_t *) shortPtr)[0];
  ((uint8_t *) shortPtr)[0] = temp;

  return *shortPtr;
}


static void
SwapDoubleWord(float *src)
{
  uint32_t *longPtr = (uint32_t *) src;
  uint16_t temp = *((uint16_t *) (((uint8_t *) longPtr) + 2));

  ((uint8_t *) longPtr)[3] = ((uint8_t *) longPtr)[0];
  ((uint8_t *) longPtr)[2] = ((uint8_t *) longPtr)[1];
  ((uint8_t *) longPtr)[1] = ((uint8_t *) &temp)[0];
  ((uint8_t *) longPtr)[0] = ((uint8_t *) &temp)[1];
}

int ReadResampleConfigFile(string resampleConfigFilename, struct ResampleParams *resampleParams)
{
  ifstream configFile (resampleConfigFilename.c_str());
  std::string line;
  double val; 
  std::string identifier;

  if (configFile.is_open()){ 
   
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline0; 
    sline0<<line;
    sline0 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->imageType=val;
 
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline; 
    sline<<line;
    sline >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->numPyrLevels=val;

    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline1; 
    sline1<<line;
    sline1 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->pyrResampleFactor=val;
  
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline2; 
    sline2<<line;
    sline2 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->imageROILeft=val;
 
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline3; 
    sline3<<line;
    sline3 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->imageROITop= val;
  
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline4; 
    sline4<<line;
    sline4 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->imageROIRight= val;
  
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline5; 
    sline5<<line;
    sline5 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->imageROIBottom=val;
 
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline6; 
    sline6<<line;
    sline6 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->tileSize=val;

    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline7; 
    sline7<<line;
    sline7 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->tilePaddingTop= val;
  
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline8; 
    sline8<<line;
    sline8 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->tilePaddingLeft= val;
  

    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline9; 
    sline9<<line;
    sline9 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->tilePaddingRight=val;
 
    getline (configFile,line);
    //cout<<line<<endl;
    stringstream sline10; 
    sline10<<line;
    sline10 >> identifier >> val;
    //cout<<val<<endl;
    resampleParams->tilePaddingBottom=val;
    configFile.close();
    return 1;
  }
  else{
    cout << "Unable to open settings file"; 
    return 0;
  }
  
}

void PrintResampleParams(struct ResampleParams *resampleParams)
{
  cout<<"IMAGE_TYPE = "<<resampleParams->imageType<<endl;
  cout<<"NUM_PYRAMID_LEVELS = "<<resampleParams->numPyrLevels<<endl;
  cout<<"PYRAMID_RESAMPLE_FACTOR = "<<resampleParams->pyrResampleFactor<<endl;
  cout<<"IMAGE_ROI_LEFT = "<<resampleParams->imageROILeft<<endl; 
  cout<<"IMAGE_ROI_TOP = "<<resampleParams->imageROITop<<endl;
  cout<<"IMAGE_ROI_RIGHT = "<<resampleParams->imageROIRight<<endl;
  cout<<"IMAGE_ROI_BOTTOM = "<<resampleParams->imageROIBottom<<endl;
  cout<<"TILE_SIZE = "<<resampleParams->tileSize<<endl;
  cout<<"TILE_PADDING_TOP = "<<resampleParams->tilePaddingTop<<endl;
  cout<<"TILE_PADDING_LEFT = "<<resampleParams->tilePaddingLeft<<endl;
  cout<<"TILE_PADDING_RIGHT = "<<resampleParams->tilePaddingRight<<endl;
  cout<<"TILE_PADDING_BOTTOM = "<<resampleParams->tilePaddingBottom<<endl;
}

//returns a binary decomposition from LSB to MSB!!
void BinaryDecomposition(int number, int numBits, vector<int> &binDecomposition)
{
 
  for (int i = numBits-1; i >=0; i--){
       int denominator = floor(pow(2.0,i));
       int divider = number/denominator;
       number = number-divider*denominator;
       binDecomposition[i]=divider;
  }
}

vector<int> QuadDecomposition(int h_number, int v_number, int numBits)
{
  vector<int>quadDecomposition;
  quadDecomposition.resize(numBits);
  vector<int> horBinDecomposition;
  vector<int> verBinDecomposition;
  horBinDecomposition.resize(numBits);
  verBinDecomposition.resize(numBits);
  
  BinaryDecomposition(h_number, numBits, horBinDecomposition);
  BinaryDecomposition(v_number, numBits, verBinDecomposition);
  
  for (int i=0; i < numBits; i++){
    quadDecomposition[i] = horBinDecomposition[i]+2*verBinDecomposition[i];
  }

  horBinDecomposition.clear();
  verBinDecomposition.clear();

  return quadDecomposition;
}

void CreateDirectoryStructure(string root, int numPyrLevels)
{
 
    float length = pow(2.0, numPyrLevels);
    for (int h = 0; h < length; h++){
      for (int v = 0; v < length; v++){

        vector<int> quadDecomposition;
	quadDecomposition = QuadDecomposition(h, v, numPyrLevels);
        stringstream commandss; 
        commandss<<string("mkdir ")<<root;
        for (int i = 0; i < quadDecomposition.size(); i++){
	  commandss<<string("/")<<(quadDecomposition[i]+1);
	}
        //cout<<commandss.str()<<endl;
        system(commandss.str().c_str());
        quadDecomposition.clear();
      }
    }
 
}

struct TileParams MakeTilingParams(struct ResampleParams resampleParams, int fullImgWidth, int fullImgHeight)
{
  
  struct TileParams tileParams;                                                 
  int imgWidth = fullImgWidth - (resampleParams.imageROILeft + resampleParams.imageROIRight);
  int imgHeight = fullImgHeight -(resampleParams.imageROITop + resampleParams.imageROIBottom); 
  
  if (resampleParams.numPyrLevels == 0){//noPyramid
  }

  if (resampleParams.numPyrLevels < 0){//automatic pyramid generation to match the size of one tile
    float heightRatio = imgHeight/(float)resampleParams.tileSize;
    float widthRatio = imgWidth/(float)resampleParams.tileSize;
    cout<<"heightRatio = "<<heightRatio<<", widthRatio="<<widthRatio<<endl;
    float sizeRatio = heightRatio;
    if (widthRatio>sizeRatio){
      sizeRatio = widthRatio;
    }  
    cout<<"sizeratio = "<<sizeRatio<<endl;

    float maxNumPyrLevels = ceil(log2(sizeRatio));
    cout<<"maxNumPyrLevels="<<maxNumPyrLevels<<endl;

    cout<<"imgWidth="<<imgWidth<<", imgHeight="<<imgHeight<<endl;
    int padImgWidth = resampleParams.tileSize*pow(2, maxNumPyrLevels);    
    int padImgHeight = resampleParams.tileSize*pow(2, maxNumPyrLevels);
    cout<<"padImgWidth="<<padImgWidth<<", padImgHeight="<<padImgHeight<<endl;

    int horOffset = (padImgWidth-imgWidth)/2.0;
    int verOffset = (padImgHeight-imgHeight)/2.0;
    cout<<"horOffset="<<horOffset<<", verOffset="<<verOffset<<endl;
    tileParams.verOffset = verOffset;
    tileParams.horOffset = horOffset;
    tileParams.padImgWidth = padImgWidth;
    tileParams.padImgHeight = padImgHeight;
    tileParams.numPyrLevels = maxNumPyrLevels;
  }
  if (resampleParams.numPyrLevels > 0){
  }
  return tileParams;
}

#if 0
void MakeDRGPyramid(string inputFilename, struct ResampleParams resampleParams, string outputDirname)
{
   DiskImageView<PixelGray<uint8> >  initImg(inputFilename);
    GeoReference initGeo;
    read_georeference(initGeo, inputFilename);
      
    InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<PixelGray<uint8> >, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> interpInitImg
	= interpolate(edge_extend(initImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());
   
    struct TileParams tileParams = MakeTilingParams(resampleParams, initImg.cols(), initImg.rows());
   
    int padImgWidth = tileParams.padImgWidth;
    int padImgHeight = tileParams.padImgHeight;
    cout<<"numPyrLevels="<<tileParams.numPyrLevels<<endl;
    cout<<"origImgPaddingWidth="<<padImgWidth<<", origImgPaddingHeight="<<padImgHeight<<endl;

    for (int p = 0; p < tileParams.numPyrLevels+1; p++){
      CreateDirectoryStructure(outputDirname, p);
    }

    //determine the leftTop point of the original image after padding
    Vector2 leftTopPixelPadding;
    leftTopPixelPadding(0) = -tileParams.horOffset+resampleParams.imageROIRight;
    leftTopPixelPadding(1) = -tileParams.verOffset+resampleParams.imageROITop;
    Vector2 leftTopPointPadding = initGeo.pixel_to_point(leftTopPixelPadding);
    Vector2 leftTopPixel;
    leftTopPixel(0)=0;
    leftTopPixel(1)=0;
    Vector2 leftTopPoint = initGeo.pixel_to_point(leftTopPixel);
    cout<<"leftTopPoint"<<leftTopPoint<<endl;
    cout<<"leftTopPointPadding"<<leftTopPointPadding<<endl;

    Matrix<double> init_H;
    init_H = initGeo.transform();
    cout<<"init_H="<<init_H<<endl;

    float resampleFactor = 1;
  
    for (int p = 0; p <= tileParams.numPyrLevels; p++){

      stringstream ss;
      ss<<p;
           
      //ImageView<float> outImg(tileSize+1, tileSize+1);
      //ImageView<PixelGray<uint8> > outImg(resampleParams.tileSize+3, resampleParams.tileSize+3);
      ImageView<float> outImg(resampleParams.tileSize + resampleParams.tilePaddingLeft + resampleParams.tilePaddingRight, 
                              resampleParams.tileSize + resampleParams.tilePaddingTop + resampleParams.tilePaddingBottom );
      //create the outputGeo - START
      GeoReference outputGeo = initGeo;
            
      //num tiles in resampled image
      int numHorTiles = (padImgWidth*resampleFactor)/resampleParams.tileSize;
      int numVerTiles = (padImgHeight*resampleFactor)/resampleParams.tileSize;
      cout<<"pyrLevel="<<p<<", numHorTiles="<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
  
      //create the georeference for the resampled and padded image
      Matrix<double> resampled_H;
      resampled_H = initGeo.transform();
      resampled_H(0,2) = leftTopPointPadding(0);
      resampled_H(1,2) = leftTopPointPadding(1);
      resampled_H(0,0) = init_H(0,0)/resampleFactor;
      resampled_H(1,1) = init_H(1,1)/resampleFactor;
      cout<<"resample_H: "<<resampled_H<<endl;
      GeoReference resampledGeo = initGeo;
      resampledGeo.set_transform(resampled_H);
      //#if 0
      //scanning the resampled image
      for (int k = 0; k <numHorTiles; k++){
	for (int l = 0; l <numVerTiles; l++){
           
          //topleft pixel position in this tile of the resampled image
          Vector2 leftTop;
          leftTop(0) = k*resampleParams.tileSize;
	  leftTop(1) = l*resampleParams.tileSize;
     
          //coordinates in the resampled padded image
	  Vector2 point = resampledGeo.pixel_to_point(leftTop);
	 
          //create the outputGeo - START
          Matrix<double> tile_H;
          tile_H = resampledGeo.transform();
	  
          //lon = H(0,0)*i + 0*j + H(0,2)
	  //lat = 0*i + H(1,1)*j + H(1,2)
          tile_H(0,2) = point(0);
	  tile_H(1,2) = point(1);
	  tile_H(0,0) = resampled_H(0,0);
	  tile_H(1,1) = resampled_H(1,1);
          cout<<"tile_H: "<<tile_H<<endl;
	  outputGeo.set_transform(tile_H);
	  //create the outputGeo - END

          Vector2 leftTopPaddingImage, rightBottomPaddingImage;
          leftTopPaddingImage(0) = leftTop(0) - resampleParams.tilePaddingLeft;
	  leftTopPaddingImage(1) = leftTop(1) - resampleParams.tilePaddingTop;
          rightBottomPaddingImage(0) = leftTopPaddingImage(0) + resampleParams.tileSize + resampleParams.tilePaddingRight;
	  rightBottomPaddingImage(1) = leftTopPaddingImage(1) + resampleParams.tileSize + resampleParams.tilePaddingBottom;

	  //int leftTopImage = leftTop(0) + resampleParams.tileSize+resampleParams.tilePaddingRight;
	  //int rightBottomImage = leftTop(1) + resampleParams.tileSize+resampleParams.tilePaddingBottom;
          //int endLeftImage = leftTop(0) + resampleParams.tileSize+resampleParams.tilePaddingRight;
	  //int endBottomImage = leftTop(1) + resampleParams.tileSize+resampleParams.tilePaddingBottom;
          
          //create a new tile of the resampled image - START
	  /*
	  for (int i = leftTop(0); i < leftTop(0) + resampleParams.tileSize+resampleParams.tilePaddingRight; i++){
	    for (int j = leftTop(1); j < leftTop(1) + resampleParams.tileSize+resampleParams.tilePaddingBottom; j++){
	  */
	  for (int i = leftTopPaddingImage(0); i < rightBottomPaddingImage(0); i++){
	    for (int j = leftTopPaddingImage(1); j < rightBottomPaddingImage(1); j++){    
              
              Vector2 leftTopTile;
              leftTopTile(0) = i-leftTopPaddingImage(0);
              leftTopTile(1) = j-leftTopPaddingImage(1);
              //int leftTopTile = i-leftTopPaddingImage(0);
              //int startLeftTile = j-leftTopPaddingImage(1);

              //pixel in the resampled image
              Vector2 resampledPix;
	      resampledPix(0) = i;
	      resampledPix(1) = j;
              Vector2 lonlat = resampledGeo.pixel_to_lonlat(resampledPix);
	      //pixel in the original image
              Vector2 initPix = initGeo.lonlat_to_pixel(lonlat);	      
	      /*
              if ((initPix(0) >=0 ) && (initPix(0) < initImg.cols()) && (initPix(1)>=0) && (initPix(1) < initImg.rows())){
		outImg.impl()(i-leftTop(0), j-leftTop(1)) = interpInitImg.impl()(initPix(0), initPix(1));
	      }
              else{
		outImg.impl()(i-leftTop(0), j-leftTop(1)) = -3.4028226550889e+38;
              }
	      */
              if ((initPix(0) >=0 ) && (initPix(0) < initImg.cols()) && (initPix(1)>=0) && (initPix(1) < initImg.rows())){
		outImg.impl()(leftTopTile(0), leftTopTile(1)) = interpInitImg.impl()(initPix(0), initPix(1));
	      }
              else{
		outImg.impl()(leftTopTile(0), leftTopTile(1)) = -3.4028226550889e+38;
              }
	    }
	  }      
          //create a new tile of the resampled image - END
          
          //create the out filenames - START     
	  //vector<int> quadDecomposition = QuadDecomposition(77/*k*/, 73/*l*/, 7/*tileParams.numPyrLevels*/);
          vector<int> quadDecomposition = QuadDecomposition(k, l, tileParams.numPyrLevels-p);
          
          //set in a directory format string
          stringstream currDir;
          for (int index=0;index < quadDecomposition.size(); index++){
            currDir<<string("/"); 
	    currDir<<(quadDecomposition[quadDecomposition.size()-1-index]+1);
          }
          //cout<<bilFilename<<endl;
          quadDecomposition.clear();
     
          stringstream horTileIndex;
	  horTileIndex<<k;
          stringstream verTileIndex;
	  verTileIndex<<l;

	  string tifFilename =  outputDirname+"/res_"+ss.str()+"_"+horTileIndex.str()+"_"+verTileIndex.str()+".tif";
	  cout<<"tifFilename="<<tifFilename<<endl;
          string bilFilename =  outputDirname+currDir.str()+"/0";
	  cout<<"bilFilename="<<bilFilename<<endl;
	  //create the output filename - END 
	  
	  //write the output .tif file
	  //write_georeferenced_image(tifFilename,
	  //			    outImg,
 	  //			    outputGeo, TerminalProgressCallback("photometry","Processing:"));

          
          int length = outImg.cols()*outImg.rows(); 
          FILE* bilfile = fopen (bilFilename.c_str(), "wb");
	  for (int jx = 0; jx < outImg.rows(); jx++){
             for (int ix = 0; ix < outImg.cols(); ix++){
	       unsigned char f = outImg(ix, jx);
               //cout<<"resample: "<<f<<endl;
	       fwrite(&f, sizeof(unsigned char), 1, bilfile);
	     }
	  }
	  fclose(bilfile);
	  
	  
	}
      }
      //#endif           
      resampleFactor = resampleFactor*resampleParams.pyrResampleFactor;
    }
}


void MakeDEMPyramid(string inputFilename, struct ResampleParams resampleParams, string outputDirname)
{
  
    DiskImageView<float>  initImg(inputFilename);
    GeoReference initGeo;
    read_georeference(initGeo, inputFilename);
      
    InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<float>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> interpInitImg
    	= interpolate(edge_extend(initImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

    struct TileParams tileParams = MakeTilingParams(resampleParams, initImg.cols(), initImg.rows());
   
    int padImgWidth = tileParams.padImgWidth;
    int padImgHeight = tileParams.padImgHeight;
    cout<<"numPyrLevels="<<tileParams.numPyrLevels<<endl;
    cout<<"origImgPaddingWidth="<<padImgWidth<<", origImgPaddingHeight="<<padImgHeight<<endl;

    for (int p = 0; p < tileParams.numPyrLevels+1; p++){
      CreateDirectoryStructure(outputDirname, p);
    }

    //determine the leftTop point of the original image after padding
    Vector2 leftTopPixelPadding;
    leftTopPixelPadding(0) = -tileParams.horOffset+resampleParams.imageROIRight;
    leftTopPixelPadding(1) = -tileParams.verOffset+resampleParams.imageROITop;
    Vector2 leftTopPointPadding = initGeo.pixel_to_point(leftTopPixelPadding);
    Vector2 leftTopPixel;
    leftTopPixel(0)=0;
    leftTopPixel(1)=0;
    Vector2 leftTopPoint = initGeo.pixel_to_point(leftTopPixel);
    cout<<"leftTopPoint"<<leftTopPoint<<endl;
    cout<<"leftTopPointPadding"<<leftTopPointPadding<<endl;

    Matrix<double> init_H;
    init_H = initGeo.transform();
    cout<<"init_H="<<init_H<<endl;

    float resampleFactor = 1;
  
    for (int p = 0; p <= tileParams.numPyrLevels; p++){

      stringstream ss;
      ss<<p;
           
      //ImageView<float> outImg(resampleParams.tileSize+1, resampleParams.tileSize+1);
      ImageView<float> outImg(resampleParams.tileSize + resampleParams.tilePaddingLeft + resampleParams.tilePaddingRight, 
                              resampleParams.tileSize+resampleParams.tilePaddingTop + resampleParams.tilePaddingBottom );
  
      //create the outputGeo - START
      GeoReference outputGeo = initGeo;
            
      //num tiles in resampled image
      int numHorTiles = (padImgWidth*resampleFactor)/resampleParams.tileSize;
      int numVerTiles = (padImgHeight*resampleFactor)/resampleParams.tileSize;
      cout<<"pyrLevel="<<p<<", numHorTiles="<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
  
      //create the georeference for the resampled and padded image
      Matrix<double> resampled_H;
      resampled_H = initGeo.transform();
      resampled_H(0,2) = leftTopPointPadding(0);
      resampled_H(1,2) = leftTopPointPadding(1);
      resampled_H(0,0) = init_H(0,0)/resampleFactor;
      resampled_H(1,1) = init_H(1,1)/resampleFactor;
      cout<<"resample_H: "<<resampled_H<<endl;
      GeoReference resampledGeo = initGeo;
      resampledGeo.set_transform(resampled_H);
      //#if 0
      //scanning the resampled image
      for (int k = 0; k <numHorTiles; k++){
	for (int l = 0; l <numVerTiles; l++){
           
          //topleft pixel position in this tile of the resampled image
          Vector2 leftTop;
          leftTop(0) = k*resampleParams.tileSize;
	  leftTop(1) = l*resampleParams.tileSize;
     
          //coordinates in the resampled padded image
	  Vector2 point = resampledGeo.pixel_to_point(leftTop);
	 
          //create the outputGeo - START
          Matrix<double> tile_H;
          tile_H = resampledGeo.transform();
	  
          //lon = H(0,0)*i + 0*j + H(0,2)
	  //lat = 0*i + H(1,1)*j + H(1,2)
          tile_H(0,2) = point(0);
	  tile_H(1,2) = point(1);
	  tile_H(0,0) = resampled_H(0,0);
	  tile_H(1,1) = resampled_H(1,1);
          cout<<"tile_H: "<<tile_H<<endl;
	  outputGeo.set_transform(tile_H);
	  //create the outputGeo - END

	  //#if 0
	  Vector2 leftTopPaddingImage, rightBottomPaddingImage;
          leftTopPaddingImage(0) = leftTop(0) - resampleParams.tilePaddingLeft;
	  leftTopPaddingImage(1) = leftTop(1) - resampleParams.tilePaddingTop;
          rightBottomPaddingImage(0) = leftTopPaddingImage(0) + resampleParams.tileSize + resampleParams.tilePaddingRight;
	  rightBottomPaddingImage(1) = leftTopPaddingImage(1) + resampleParams.tileSize + resampleParams.tilePaddingBottom;

	  for (int i = leftTopPaddingImage(0); i < rightBottomPaddingImage(0); i++){
	    for (int j = leftTopPaddingImage(1); j < rightBottomPaddingImage(1); j++){    
              
              Vector2 leftTopTile;
              leftTopTile(0) = i-leftTopPaddingImage(0);
              leftTopTile(1) = j-leftTopPaddingImage(1);
           
              //pixel in the resampled image
              Vector2 resampledPix;
	      resampledPix(0) = i;
	      resampledPix(1) = j;
              Vector2 lonlat = resampledGeo.pixel_to_lonlat(resampledPix);
	      //pixel in the original image
              Vector2 initPix = initGeo.lonlat_to_pixel(lonlat);	      
	  
              if ((initPix(0) >=0 ) && (initPix(0) < initImg.cols()) && (initPix(1)>=0) && (initPix(1) < initImg.rows())){
		outImg.impl()(leftTopTile(0), leftTopTile(1)) = interpInitImg.impl()(initPix(0), initPix(1));
	      }
              else{
		outImg.impl()(leftTopTile(0), leftTopTile(1)) = -3.4028226550889e+38;
              }
	    }
	  }      
          //create a new tile of the resampled image - END
	  //#endif

#if 0
	  //this is obsolete - START
	  for (int i = leftTop(0); i < leftTop(0) + resampleParams.tileSize+resampleParams.tilePaddingRight; i++){
	    for (int j = leftTop(1); j < leftTop(1) + resampleParams.tileSize+resampleParams.tilePaddingBottom; j++){
	      
              //pixel in the resampled image
              Vector2 resampledPix;
	      resampledPix(0) = i;
	      resampledPix(1) = j;
              Vector2 lonlat = resampledGeo.pixel_to_lonlat(resampledPix);
	      //pixel in the original image
              Vector2 initPix = initGeo.lonlat_to_pixel(lonlat);	      

              if ((initPix(0) >=0 ) && (initPix(0) < initImg.cols()) && (initPix(1)>=0) && (initPix(1) < initImg.rows())){
		outImg.impl()(i-leftTop(0), j-leftTop(1)) = interpInitImg.impl()(initPix(0), initPix(1));
	      }
              else{
		outImg.impl()(i-leftTop(0), j-leftTop(1)) = -3.4028226550889e+38;
              }
            
	    }
	  }     
           //this is obsolete - END 
#endif
          //create a new tile of the resampled image - END
          
          //create the out filenames - START     
	  //vector<int> quadDecomposition = QuadDecomposition(77/*k*/, 73/*l*/, 7/*tileParams.numPyrLevels*/);
          vector<int> quadDecomposition = QuadDecomposition(k, l, tileParams.numPyrLevels-p);
          
          //set in a directory format string
          stringstream currDir;
          for (int index=0;index < quadDecomposition.size(); index++){
            currDir<<string("/"); 
	    currDir<<(quadDecomposition[quadDecomposition.size()-1-index]+1);
          }
          quadDecomposition.clear();
     
          stringstream horTileIndex;
	  horTileIndex<<k;
          stringstream verTileIndex;
	  verTileIndex<<l;

	  string tifFilename =  outputDirname+"/res_"+ss.str()+"_"+horTileIndex.str()+"_"+verTileIndex.str()+".tif";
	  cout<<"tifFilename="<<tifFilename<<endl;
          string bilFilename =  outputDirname+currDir.str()+"/0";
	  cout<<"bilFilename="<<bilFilename<<endl;
	  //create the output filename - END 
	  
	  //write the output .tif file
	  /*
	  write_georeferenced_image(tifFilename,
				    outImg,
				    outputGeo, TerminalProgressCallback("photometry","Processing:"));
	  */
          //write the output .bil file
          int length = outImg.cols()*outImg.rows(); 
          FILE* bilfile = fopen (bilFilename.c_str(), "wb");
	  for (int jx = 0; jx < outImg.rows(); jx++){
             for (int ix = 0; ix < outImg.cols(); ix++){
	       float f = outImg(ix, jx);
               SwapDoubleWord(&f);
	       fwrite(&f, sizeof(float), 1, bilfile);
	     }
	  }
	  fclose(bilfile);
	}
      }
      //#endif           
      resampleFactor = resampleFactor*resampleParams.pyrResampleFactor;
    }
}
#endif

//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) {
    
    string inputFilename;
    string outputDirname;
    string configFilename;

    po::options_description general_options("Options");

    general_options.add_options()
    ("inputFilename,i", po::value<std::string>(&inputFilename)->default_value("input.tif"), "input filename.")
    ("outputDirname,o", po::value<std::string>(&outputDirname)->default_value("result"), "output dirname.")
    ("configFilename,c", po::value<std::string>(&configFilename)->default_value("resample_config.txt"), "configuration filename.")
    ("help,h", "Display this help message");
    po::options_description hidden_options("");
       
    po::options_description options("Allowed Options");
    options.add(general_options);
    
    po::positional_options_description p;
    p.add("inputFilename", -1);
    
    std::ostringstream usage;
    usage << "Description: main code for georeferenced file resample and tiling" << std::endl << std::endl;
    usage << general_options << std::endl;
    
    po::variables_map vm;
    try {
      po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
      po::notify( vm );
    } catch ( po::error &e ) {
      std::cout << "An error occured while parsing command line arguments.\n";
      std::cout << "\t" << e.what() << "\n\n";
      std::cout << usage.str();
      return 1;
    }
    
    if( vm.count("help") ) {
      std::cerr << usage.str() << std::endl;
      return 1;
    }
    
    if( vm.count("inputFilename") < 1 ) {
      std::cerr << "Error: Must specify at least one file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
    }
    
    //resampling factor: < 1 is downsampling the image, > 1 is upsampling the image   

    struct ResampleParams resampleParams;
    int error = ReadResampleConfigFile(configFilename, &resampleParams);
    PrintResampleParams(&resampleParams);
  
    
    cout<<"inputFilename="<<inputFilename<<endl;
    boost::shared_ptr<DiskImageResource> img_rsrc( new DiskImageResourceGDAL(inputFilename) );
    if (img_rsrc->has_nodata_read()){
      float noDataVal = img_rsrc->nodata_read();
      cout<<"noDataVal="<<noDataVal;
    }
    else{
      cout<<"no Data Val not found"<<endl;
    }
    /*
    if (resampleParams.imageType == 0){
      //MakeDEMPyramid(inputFilename, resampleParams, outputDirname);
      DiskImageView<float >  initImg(inputFilename);
      GeoReference initGeo;
      read_georeference(initGeo, inputFilename);
      MakePyramid(initImg, initGeo, resampleParams, outputDirname);
    }
    
    if (resampleParams.imageType == 1){     
      //MakeDRGPyramid(inputFilename, resampleParams, outputDirname);
      DiskImageView<PixelGray<uint8> >  initImg(inputFilename);
      GeoReference initGeo;
      read_georeference(initGeo, inputFilename);
      MakePyramid(initImg, initGeo, resampleParams, outputDirname);
    } 
    */

}


















