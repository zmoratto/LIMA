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

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;

float ComputePixPerDegree(GeoReference Geo, int width, int height, int useUSGS_lonlat)
{

  //cout << Geo <<"\n";
  float radius = Geo.datum().semi_major_axis();
  cout<<"radius="<<radius<<endl;
  //float meridianOffset = Geo.datum().meridian_offset();
  //cout<<"meridianOffset="<<meridianOffset<<endl;
 
  printf("width = %d, height = %d\n", width, height);
 
  Vector2 leftTopPixel(0,0);
  Vector2 leftTopLonLat = Geo.pixel_to_lonlat(leftTopPixel);

  Vector2 rightBottomPixel(width-1, height-1);
  Vector2 rightBottomLonLat = Geo.pixel_to_lonlat(rightBottomPixel);
  
  float minLon = leftTopLonLat(0);
  float minLat = leftTopLonLat(1);
  float maxLon = rightBottomLonLat(0);
  float maxLat = rightBottomLonLat(1);
  
  printf("BEFORE: minLon=%f, maxLon=%f, minLat=%f, maxLat=%f\n", 
                  minLon, maxLon, minLat, maxLat);  
 
  
  if (useUSGS_lonlat == 1){
     float usgs_2_lonlat = 180/(3.14159265*3396190);
     minLon = minLon*usgs_2_lonlat; 
     maxLon = maxLon*usgs_2_lonlat; 
     minLat = minLat*usgs_2_lonlat; 
     maxLat = maxLat*usgs_2_lonlat;
  }

  if (minLat < 0) minLat = minLat+180;
  if (minLon < 0) minLon = minLon+180;
  if (maxLat < 0) maxLat = maxLat+180;
  if (maxLon < 0) maxLon = maxLon+180;
  
  printf("AFTER: minLon=%f, maxLon=%f, minLat=%f, maxLat=%f\n",
  	         minLon, maxLon, minLat, maxLat);

  float numPixPerDegree;
  numPixPerDegree = width/fabs(maxLon-minLon);
  numPixPerDegree = height/fabs(maxLat-minLat);
  
  return numPixPerDegree;
}


int main( int argc, char *argv[] ) {
    
    string initFilename;
    string correctFilename;
    string correctedFilename;

    po::options_description general_options("Options");

    general_options.add_options()
    ("initFilename,i", po::value<std::string>(&initFilename)->default_value("init"), "initial filename.")
    ("correctFilename,c", po::value<std::string>(&correctFilename)->default_value("correct"), "correct filename.")
    ("correctedFilename,r", po::value<std::string>(&correctedFilename)->default_value("result"), "result filename.")
    ("help,h", "Display this help message");
    po::options_description hidden_options("");
       
    po::options_description options("Allowed Options");
    options.add(general_options);
    
    po::positional_options_description p;
    p.add("initFilenme", -1);
    
    std::ostringstream usage;
    usage << "Description: main code for georeference correction" << std::endl << std::endl;
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
    
    if( vm.count("initFilename") < 1 ) {
      std::cerr << "Error: Must specify at least one file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
    }
       
    DiskImageView<PixelGray<uint8> >  initFile(initFilename);
    GeoReference initGeo;
    read_georeference(initGeo, initFilename);

    float init_radius = initGeo.datum().semi_major_axis();
    cout<<"init_radius="<<init_radius<<endl;

    DiskImageView<float>   correctFile(correctFilename);
    GeoReference correctGeo;
    read_georeference(correctGeo, correctFilename);
 
 
    //create the correctedGeo - START
    GeoReference correctedGeo = correctGeo;
    float correct_radius = correctedGeo.datum().semi_major_axis();
    cout<<"correct_radius="<<correct_radius<<endl;

    Matrix<double> init_H;
    init_H = initGeo.transform();
    cout<<"init_H="<<init_H<<endl;
   
    Matrix<double> correct_H;
    correct_H = correctGeo.transform();
    
    Matrix<double> corrected_H;
    corrected_H = correctedGeo.transform();
    //lon = H(0,0)*i+0*j + H(0,2)
    //lat = 0*i+H(1,1)*j + H(1,2)
    
    //corrected_H(0,2)=init_H(0,2);
    corrected_H(0,2) = init_H(0,2) - (175.5-180)*(3.14159265*correct_radius/180);
    corrected_H(1,2)=init_H(1,2);
    corrected_H(0,0)=init_H(0,0);
    corrected_H(1,1)=init_H(1,1);
    correctedGeo.set_transform(corrected_H);
    //create the correctedGeo - END

    //find  the image boundaries with the init geo 
    Vector2 leftTopPix;
    leftTopPix(0) = 0; //x coordinate
    leftTopPix(1) = 0; //y coordinate
    
    Vector2 rightBottomPix;
    rightBottomPix(0) = initFile.cols()-1; //x coordinate
    rightBottomPix(1) = initFile.rows()-1; //y coordinate

    cout<<"------------------"<<endl;
    
    cout<<"correct_H="<<correct_H<<endl;
    cout<<"corrected_H="<<corrected_H<<endl;
    
    cout<<"------------------"<<endl;
    cout<<"Cropped Image Boundaries using the init and corrected geo respectively"<<endl;
 
    Vector2 initLeftTopLonLat = initGeo.pixel_to_lonlat(leftTopPix);
    initLeftTopLonLat(0)=180+initLeftTopLonLat(0)/(3.14159265*correct_radius/180);
    initLeftTopLonLat(1)=initLeftTopLonLat(1)/(3.14159265*correct_radius/180);
    cout<<"initLeftTopLonLat="<<initLeftTopLonLat<<endl;

    Vector2 correctedLeftTopLonLat = correctedGeo.pixel_to_lonlat(leftTopPix);
    cout<<"correctedLeftTopLonLat="<<correctedLeftTopLonLat<<endl;
    
    Vector2 initRightBottomLonLat = initGeo.pixel_to_lonlat(rightBottomPix);
    initRightBottomLonLat(0)=180+initRightBottomLonLat(0)/(3.14159265*correct_radius/180);
    initRightBottomLonLat(1)=initRightBottomLonLat(1)/(3.14159265*correct_radius/180);
    cout<<"initRightBottomLonLat="<<initRightBottomLonLat<<endl;

    Vector2 correctedRightBottomLonLat = correctedGeo.pixel_to_lonlat(rightBottomPix);
    cout<<"correctedRightBottomLonLat="<<correctedRightBottomLonLat<<endl;
    
    cout<<"------------------"<<endl;
    
    Vector2 centerLonLat;
    centerLonLat(0) = 175.526;
    centerLonLat(1) = -14.6064;
    cout <<"centerLonLat="<<centerLonLat<<endl; 
    Vector2 centerPix = correctedGeo.lonlat_to_pixel(centerLonLat);
    cout <<"centerPix="<<centerPix<<endl; 

    //write the corrected file
    write_georeferenced_image(correctedFilename,
                              initFile,
                              correctedGeo, TerminalProgressCallback("photometry","Processing:"));
    

}


















