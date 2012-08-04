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
#include "string_util.h"
#include "geotif_resample.h"

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;

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



//utility to upsample/downsample and tile geotiff file 
int main( int argc, char *argv[] ) {
    
    string inputFilename;
    string outputDirname;
    string configFilename;

    int horTile, verTile, numPyrLevels;

    po::options_description general_options("Options");

    general_options.add_options()
    ("inputFilename,i", po::value<std::string>(&inputFilename)->default_value("input.tif"), "input filename.")
    ("outputDirname,o", po::value<std::string>(&outputDirname)->default_value("result"), "output dirname.")
    ("horizontalTile,x", po::value<int>(&horTile)->default_value(0), "horizontal tile.")
    ("verticalTile,y", po::value<int>(&verTile)->default_value(0), "vertical tile.")
    ("numPyramidLevels,n", po::value<int>(&numPyrLevels)->default_value(0), "number of pyramid levels.")
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
      cout<<"noDataVal not found"<<endl;
    }
    

   string inputFilenameNoPath = GetFilenameNoPath(inputFilename);  
   std::vector<std::string> strs;
   strs = FindAndSplit(inputFilenameNoPath, string("_"));
   istringstream (strs[1])>>horTile; 
   istringstream (strs[2])>>verTile;
   //cout<<"horTile="<<horTile<<", verTile="<<verTile<<endl; 

    vector<int> quadDecomposition;
    quadDecomposition = QuadDecomposition(horTile, verTile, numPyrLevels);
    for (int i = numPyrLevels-1; i >=0; i--){
      quadDecomposition[i] = quadDecomposition[i] +1; 
      stringstream ss;//create a stringstream
      ss << quadDecomposition[i];
      outputDirname = outputDirname+string("/")+ss.str();
    }
    cout<<"outputDirname="<<outputDirname<<endl;
   
  
    
    if (resampleParams.imageType == 0){//DEM
      DiskImageView<float >  initImg(inputFilename);
      GeoReference initGeo;
      read_georeference(initGeo, inputFilename);
      cout<<"make pyramid"<<endl;
      MakePyramid(initImg, initGeo, resampleParams, outputDirname);
    }
    
    if (resampleParams.imageType == 1){//DRG     
      DiskImageView<PixelGray<uint8> >  initImg(inputFilename);
      GeoReference initGeo;
      read_georeference(initGeo, inputFilename);
      MakePyramid(initImg, initGeo, resampleParams, outputDirname);
    } 
    

}


















