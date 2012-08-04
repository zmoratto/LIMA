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

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;

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

static bool CPUIsBigEndian();
static uint16_t SwapWord(unsigned short src);
static void SwapDoubleWord(float *src);

//returns a binary decomposition from LSB to MSB!!
void BinaryDecomposition(int number, int numBits, vector<int> &binDecomposition);
vector<int> QuadDecomposition(int h_number, int v_number, int numBits);
void CreateDirectoryStructure(string root, int numPyrLevels);

template <class ViewT1>
void
MakePyramid(ImageViewBase<ViewT1> const& initImg, GeoReference const &initGeo,
             struct ResampleParams resampleParams, string outputDirname)
{
    InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<typename ViewT1::pixel_type>, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> interpInitImg
    = interpolate(edge_extend(initImg.impl(),ConstantEdgeExtension()), BilinearInterpolation());

    struct TileParams tileParams = MakeTilingParams(resampleParams, initImg.impl().cols(), initImg.impl().rows());
   
    /*
    //determine the output dirname
    vector<int> quadDecomposition;
    int h = 76;
    int v = 73;
    int numPyrLevels = 7;
    quadDecomposition = QuadDecomposition(h, v, numPyrLevels);
    for (int i = numPyrLevels-1; i >=0; i--){
       quadDecomposition[i] = quadDecomposition[i] +1; 
       cout<<quadDecomposition[i]<<endl;
    }
    */
    
    int padImgWidth = tileParams.padImgWidth;
    int padImgHeight = tileParams.padImgHeight;
    //cout<<"numPyrLevels="<<tileParams.numPyrLevels<<endl;
    //cout<<"origImgPaddingWidth="<<padImgWidth<<", origImgPaddingHeight="<<padImgHeight<<endl;

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
    //cout<<"leftTopPoint"<<leftTopPoint<<endl;
    //cout<<"leftTopPointPadding"<<leftTopPointPadding<<endl;

    Matrix<double> init_H;
    init_H = initGeo.transform();
    cout<<"init_H="<<init_H<<endl;

    float resampleFactor = 1;
  
    for (int p = 0; p <= tileParams.numPyrLevels; p++){

      stringstream ss;
      ss<<p;
     
      
      ImageView<typename ViewT1::pixel_type> outImg(resampleParams.tileSize + resampleParams.tilePaddingLeft + resampleParams.tilePaddingRight, 
						    resampleParams.tileSize + resampleParams.tilePaddingTop + resampleParams.tilePaddingBottom );
         
      
      //create the outputGeo - START
      GeoReference outputGeo = initGeo;
            
      //num tiles in resampled image
      int numHorTiles = (padImgWidth*resampleFactor)/resampleParams.tileSize;
      int numVerTiles = (padImgHeight*resampleFactor)/resampleParams.tileSize;
      //cout<<"pyrLevel="<<p<<", numHorTiles="<<numHorTiles<<", numVerTiles="<<numVerTiles<<endl;
  
      //create the georeference for the resampled and padded image
      Matrix<double> resampled_H;
      resampled_H = initGeo.transform();
      resampled_H(0,2) = leftTopPointPadding(0);
      resampled_H(1,2) = leftTopPointPadding(1);
      resampled_H(0,0) = init_H(0,0)/resampleFactor;
      resampled_H(1,1) = init_H(1,1)/resampleFactor;
      //cout<<"resample_H: "<<resampled_H<<endl;
      GeoReference resampledGeo = initGeo;
      resampledGeo.set_transform(resampled_H);
      
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
          //cout<<"tile_H: "<<tile_H<<endl;
	  outputGeo.set_transform(tile_H);
	  //create the outputGeo - END

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
	   
              if ((initPix(0) >=0 ) && (initPix(0) < initImg.impl().cols()) && (initPix(1)>=0) && (initPix(1) < initImg.impl().rows())){
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

          //int num_channels = PixelNumChannels<typename ViewT1::pixel_type>::value;
	  //cout << "num_channels = "<< num_channels << "\n";

                  
          int length = outImg.cols()*outImg.rows(); 
          FILE* bilfile = fopen (bilFilename.c_str(), "wb");
          if (bilfile == NULL){
	    cout<<bilFilename<<" not found"<<endl;
          }

	  if (resampleParams.imageType == 0){//DEM
	    for (int jx = 0; jx < outImg.rows(); jx++){
	      for (int ix = 0; ix < outImg.cols(); ix++){
                float f = outImg(ix, jx);
                SwapDoubleWord(&f);
		fwrite(&f, sizeof(float), 1, bilfile);
	      }
	    }
	  }

          if (resampleParams.imageType == 1){//DRG
	    for (int jx = 0; jx < outImg.rows(); jx++){
	      for (int ix = 0; ix < outImg.cols(); ix++){
		unsigned char f = outImg(ix, jx);
		fwrite(&f, sizeof(unsigned char), 1, bilfile);
	      }
	    }
	  }
       
	  fclose(bilfile);
      
	}
      }
     
      resampleFactor = resampleFactor*resampleParams.pyrResampleFactor;
    }

}



















