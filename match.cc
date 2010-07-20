// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Photometry.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
//#include <cv.h>
//#include <highgui.h>

#include "io.h"
#include "coregister.h"

float ComputeScaleFactor(vector<Vector3> allImgPts, vector<float> reflectance)
{
  float nominator =0.0;
  int numValidPts = 0;
  float scaleFactor = 1;

  for(int m = 0; m < reflectance.size(); m ++){
    if ((reflectance[m]!= -1) && (reflectance[m] !=0.0)){
      nominator += allImgPts[m][0]/reflectance[m];
      numValidPts += 1;
    }
  }

  if (numValidPts != 0){ 
    scaleFactor = nominator/numValidPts;
  }

  return scaleFactor;
}

vector<float> ComputeSyntImgPts(float scaleFactor, vector<float> reflectance)
{
  vector<float>synthImg;
  synthImg.resize(reflectance.size());
  for(int m = 0; m < reflectance.size(); m ++){
    if (reflectance[m] == -1){
      synthImg[m] = -1;
    }
    else{
      synthImg[m] = reflectance[m]*scaleFactor;
    }
  }

  return synthImg;
}

float ComputeMatchingError(vector<float> reflectancePts, vector<float>imgPts)
{
  float error = 0.0;
  for (int i = 0; i < reflectancePts.size(); i++){
    if (reflectancePts[i]!=-1){
      error = error + fabs(reflectancePts[i]-imgPts[i]);
    }
  }
  return error;
}


vector<Vector3>  GetAllPtsFromImage(vector<vector<LOLAShot > > trackPts,string DRGFilename)
{ 
  //vector<Vector3> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename)
  DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  vector<Vector3> allImgPts;

  int num_allImgPts = 0;
  for(int k = 0; k<trackPts.size();k++)
  {
    num_allImgPts += trackPts[k].size();
  }

  allImgPts.resize(num_allImgPts);  
  vector<pointCloud> ptHere;

  ImageViewRef<PixelMask<PixelGray<uint8> > >  interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      ptHere = trackPts[k][i].LOLAPt;
      pointCloud centerPt  = ::GetPointFromIndex( ptHere, 3);
      pointCloud topPt     = GetPointFromIndex( ptHere, 2);
      pointCloud leftPt    = GetPointFromIndex( ptHere, 1);

      if((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){

        float lon = centerPt.coords[0];
        float lat = centerPt.coords[1];
        float rad = centerPt.coords[2];

        Vector2 DEM_lonlat(lon, lat);
        Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

        int x = (int)DRG_pix[0];
        int y = (int)DRG_pix[1];

        PixelMask<PixelGray<uint8> > DRGVal = interpDRG(x, y);

        //insert data
        allImgPts[i][0] = (float) DRGVal;
        allImgPts[i][1] = x;
        allImgPts[i][2] = y;
      }else {
        //write -1 to designate and invalid point
        allImgPts[i][0] = -1;
        allImgPts[i][1] = -1;
        allImgPts[i][2] = -1;
      }
    }
  }
  return allImgPts;
}
/*
float ComputeReflectance(vector<pointCloud> LOLAPts, ModelParams modelParams, GlobalParams globalParams)
{
  pointCloud centerPt = GetPointFromIndex(LOLAPts, 3);
  pointCloud topPt = GetPointFromIndex(LOLAPts, 2);
  pointCloud leftPt = GetPointFromIndex(LOLAPts, 1);

  if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){
    Datum moon;
    moon.set_well_known_datum("D_MOON");

    centerPt.coords[2] = (centerPt.coords[2]-1737.4)*1000;
    topPt.coords[2] = (topPt.coords[2]-1737.4)*1000;
    leftPt.coords[2] = (leftPt.coords[2]-1737.4)*1000;
    printf("c = %f, t = %f, l = %f\n", centerPt.coords[2],topPt.coords[2],leftPt.coords[2]);

    Vector3 xyz = moon.geodetic_to_cartesian(centerPt.coords);
    Vector3 xyzTop = moon.geodetic_to_cartesian(topPt.coords);
    Vector3 xyzLeft = moon.geodetic_to_cartesian(leftPt.coords);
    printf("xyz = %f %f %f\n", xyz(0), xyz(1), xyz(2));
    Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);
    //printf("normal = %f %f %f\n", normal(0), normal(1), normal(2));
    float reflectance = ComputeReflectance(normal, xyz, modelParams, globalParams);

    return reflectance;
  }
  else{
    return -1;
  }
}



vector<float> ComputeAllReflectance(vector< vector<LOLAShot> > trackPts, ModelParams modelParams, GlobalParams globalParams)
{
  vector<float> reflectance;
  int num_rflc = 0;
  for(int k = 0; k < trackPts.size(); k++) 
  {
    num_rflc += trackPts[k].size() ;
  }
  
  reflectance.resize(num_rflc);
  vector<LOLAShot> LolaPts;
  for(int k = 0; k < trackPts.size();k++){
    for (int m = 0; m < trackPts[k].size();m++){
      LolaPts = trackPts[k][m];
      reflectance[m] = ComputeReflectance(LolaPts, modelParams, globalParams);
     // printf("ref = %f\n", reflectance[m]);
    }
  }
  return reflectance;
}
*/
/*
float ComputeReflectance(vector<pointCloud> LOLAPts, ModelParams modelParams, GlobalParams globalParams)
{
  pointCloud centerPt = GetPointFromIndex(LOLAPts, 3);
  pointCloud topPt = GetPointFromIndex(LOLAPts, 2);
  pointCloud leftPt = GetPointFromIndex(LOLAPts, 1);

  if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){
    Datum moon;
    moon.set_well_known_datum("D_MOON");

    centerPt.coords[2] = (centerPt.coords[2]-1737.4)*1000;
    topPt.coords[2] = (topPt.coords[2]-1737.4)*1000;
    leftPt.coords[2] = (leftPt.coords[2]-1737.4)*1000;
    printf("c = %f, t = %f, l = %f\n", centerPt.coords[2],topPt.coords[2],leftPt.coords[2]);

    Vector3 xyz = moon.geodetic_to_cartesian(centerPt.coords);
    Vector3 xyzTop = moon.geodetic_to_cartesian(topPt.coords);
    Vector3 xyzLeft = moon.geodetic_to_cartesian(leftPt.coords);
    printf("xyz = %f %f %f\n", xyz(0), xyz(1), xyz(2));
    Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);
    //printf("normal = %f %f %f\n", normal(0), normal(1), normal(2));
    float reflectance = ComputeReflectance(normal, xyz, modelParams, globalParams);

    return reflectance;
  }
  else{
    return -1;
  }
}
vector<float> ComputeTrackReflectance(vector<LOLAShot> trackPts, ModelParams modelParams, GlobalParams globalParams)
{
  vector<float> reflectance;
  reflectance.resize(trackPts.size());
  for (int m = 0; m < trackPts.size();m++){
    reflectance[m] = ComputeReflectance(trackPts[m].LOLAPt, modelParams, globalParams);
    printf("ref = %f\n", reflectance[m]);
  }
  return reflectance;
}

*/

vector<float> ComputeAllReflectance( vector< vector<LOLAShot> >  allTracks, ModelParams modelParams, GlobalParams globalParams)
{
  vector<float> reflectance;
  vector<float> single_track_ref;
  int num_reflc = 0;

  for(int k = 0; k < allTracks.size(); k++ )
  {
    num_reflc += allTracks[k].size(); 
  }
  cout << "num_reflc: " << num_reflc << endl;

  //reflectance.resize();
  cout << "reflectance.size() should be: " << reflectance.size() << endl; 
  for( int k = 0; k < allTracks.size(); k++)
  {    
    // call computeReflectance
    single_track_ref = ComputeTrackReflectance( allTracks[k], modelParams, globalParams);

    // read into vector
    for(int i = 0; i < single_track_ref.size(); i++)
    {
      reflectance.push_back(single_track_ref[i]);
    }

  }
  //read out counters - keep the party rolling
  printf("reflectance should be: %d, reflectance size is: %d\n", num_reflc, reflectance.size() );

  //return data
  return reflectance;
}



Vector<float,6> UpdateMatchingParams(vector<vector<LOLAShot> > trackPts, string DRGFilename, ModelParams modelParams,GlobalParams globalParams)
{
  DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);

  //get the true image points
  cout << "UMP calling: GetTrackPtsFromImage..." << endl;
  vector<Vector3> imgPts;
  imgPts = GetAllPtsFromImage( trackPts, DRGFilename );

  //compute the synthetic image values
  cout << "UMP: ComputeTrackReflectance..." << endl;
  vector<float> reflectance = ComputeAllReflectance( trackPts, modelParams, globalParams);
  
  cout << "UMP: ComputeScaleFactor..." << endl;
  float scaleFactor = ComputeScaleFactor(imgPts, reflectance);
  //isn't this a bit self-serving?  The scale factor should be computed a priori or conditioned from a learned distribution.  Here, we confuse training & test data.

  cout << "UMP: ComputeSynthImgPts..." << endl;
  vector<float> synthImg = ComputeSyntImgPts(scaleFactor, reflectance);


  cout << "UMP: derivative_filter..."<< endl; 
  //ImageView<float> x_deriv = derivative_filter(DRG, 1, 0);
  //ImageView<float> y_deriv = derivative_filter(DRG, 0, 1);

  ImageView<float> x_deriv;
  ImageView<float> y_deriv;

  if( DRGFilename == "../../data/Apollo15-DRG/1134_1135-DRG.tif" )
  {
    // load dirivatives that have already been calculated
    string x_dirv_load_name = "../results/Apennine_escarpment_x-dir.tiff";
    string y_dirv_load_name = "../results/Apennine_escarpment_y-dir.tiff" ;

    read_image( x_deriv, x_dirv_load_name);
    read_image( y_deriv, y_dirv_load_name);
    printf("UMP: loaded x_deriv & y_deriv...\n");
  }else{
    //calculate derivative image & save under the correct name
    string x_dirv_name = "../results/Apennine_escarpment_x-dir.tiff";
    string y_dirv_name = "../results/Apennine_escarpment_y-dir.tiff";   
    x_deriv = derivative_filter( DRG, 1, 0);
    y_deriv = derivative_filter( DRG, 0, 1);
      
    write_image( x_dirv_name, x_deriv, TerminalProgressCallback("vw", "saving x_deriv: "));
    write_image( y_dirv_name, y_deriv, TerminalProgressCallback("vw", "saving y_deriv: "));
    
  }


  cout << "UMP: calculate min & max pixel locations..."<< endl;
  int i_min = 0;
  int i_max = x_deriv.rows();
  int j_min = 0;
  int j_max = y_deriv.cols();
  printf("x range from: %f -> %f \ny range from: %f -> %f \n",i_min,i_max,j_min,j_max);

  //the following is for affine (NOT perspective) transforms
  Vector<float,6> d;//defines the affine transform
  d[0] = 1.0; d[1] = 0.0; d[2] = 0.0;
  d[3] = 0.0; d[4] = 1.0; d[5] = 0.0;
  printf("d = [ %d, %d, %d, %d %d %d]\n",d[0],d[1],d[2],d[3],d[4],d[5]);
  Matrix<float,6,6> rhs;
  Vector<float,6> lhs;

  int ii, jj;
  float x_base, y_base;
  float xx,yy;
  int row_max, col_max;
  
  //get row_max, col_max
  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
 /*
  //x_base and y_base are the initial 2D positions of the image points
  x_base = imgPts[0][1];
  y_base = imgPts[0][2];
  //how does the following make sense - ii & jj are not initialized?/!
  float xx = x_base + d[0] * ii + d[1] * jj + d[2];
  float yy = y_base + d[3] * ii + d[4] * jj + d[5];
   */
  cout << "UMP: interpolate ..." << endl ;
  InterpolationView<EdgeExtensionView<DiskImageView<PixelMask<PixelGray<uint8> > > , ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
    interpolate(DRG, BilinearInterpolation(), ZeroEdgeExtension());
  int iter=0;
  cout << "line 162, central iteration UpdateModelParams" << endl; 
  while(iter < 2) //gradient descent => optimal transform
  {
    int num_valid = 0;
    for (int i = 0; i < imgPts.size(); i++)
    {
      if( (synthImg[i]!= -1) && (synthImg[i]!=0) )
      {
        num_valid += 1;
        // lola pixel coordinates
        x_base = imgPts[i][1];
        y_base = imgPts[i][2];

        // image coordinates under the current transform
        ii = (int) floor(d[0]*x_base + d[1]*y_base + d[2]);
        jj = (int) floor(d[3]*x_base + d[4]*y_base + d[5]);
       
        if(i<50){
          printf("LOLA pnt: %d of intensity %f, at (%d,%d) transformed to (%d,%d)\n",i,imgPts[i][0],x_base,y_base,ii,jj);
        }

        // check (ii,jj) are inside the image!
        if( ( ii> 0) && ( ii<= x_deriv.rows()) && ( jj> 0) && ( jj<= y_deriv.cols())){

        //initialize constants
        float I_x_sqr, I_x_I_y, I_y_sqr; 
        float I_e_val = imgPts[i][0]-synthImg[i];
        float I_y_val, I_x_val;

        //calculate numerical dirivatives (ii,jj)... 
        if( i <100)
        {
        //  printf("i = %d: (%d,%d)\n",i,ii,jj);
        }

        I_x_val = x_deriv(ii,jj); 
        I_y_val = y_deriv(ii,jj); 
        I_x_I_y = I_x_val*I_y_val;        
        I_x_sqr = I_x_sqr*I_x_sqr;
        I_y_sqr = I_y_sqr*I_y_sqr;

        // Left hand side
        lhs(0) += ii * I_x_val * I_e_val;
        lhs(1) += jj * I_x_val * I_e_val;
        lhs(2) +=      I_x_val * I_e_val;
        lhs(3) += ii * I_y_val * I_e_val;
        lhs(4) += jj * I_y_val * I_e_val;
        lhs(5) +=      I_y_val * I_e_val;

        // Right Hand Side UL
        rhs(0,0) += ii*ii * I_x_sqr;
        rhs(0,1) += ii*jj * I_x_sqr;
        rhs(0,2) += ii    * I_x_sqr;
        rhs(1,1) += jj*jj * I_x_sqr;
        rhs(1,2) += jj    * I_x_sqr;
        rhs(2,2) +=         I_x_sqr;

        // Right Hand Side UR
        rhs(0,3) += ii*ii * I_x_I_y;
        rhs(0,4) += ii*jj * I_x_I_y;
        rhs(0,5) += ii    * I_x_I_y;
        rhs(1,4) += jj*jj * I_x_I_y;
        rhs(1,5) += jj    * I_x_I_y;
        rhs(2,5) +=         I_x_I_y;

        // Right Hand Side LR
        rhs(3,3) += ii*ii * I_y_sqr;
        rhs(3,4) += ii*jj * I_y_sqr;
        rhs(3,5) += ii    * I_y_sqr;
        rhs(4,4) += jj*jj * I_y_sqr;
        rhs(4,5) += jj    * I_y_sqr;
        rhs(5,5) +=         I_y_sqr;
        }// end of if statement: inside image
      }// end of if statement: valid reflectance  
    }// end of for loop over all data  points
    cout << "num_valid: "<< num_valid << endl;

    // Fill in symmetric entries
    rhs(1,0) = rhs(0,1);
    rhs(2,0) = rhs(0,2);
    rhs(2,1) = rhs(1,2);
    rhs(1,3) = rhs(0,4);
    rhs(2,3) = rhs(0,5);
    rhs(2,4) = rhs(1,5);
    rhs(3,0) = rhs(0,3);
    rhs(3,1) = rhs(1,3);
    rhs(3,2) = rhs(2,3);
    rhs(4,0) = rhs(0,4);
    rhs(4,1) = rhs(1,4);
    rhs(4,2) = rhs(2,4);
    rhs(4,3) = rhs(3,4);
    rhs(5,0) = rhs(0,5);
    rhs(5,1) = rhs(1,5);
    rhs(5,2) = rhs(2,5);
    rhs(5,3) = rhs(3,5);
    rhs(5,4) = rhs(4,5);
    try {
      solve_symmetric_nocopy(rhs,lhs);
    } catch (ArgumentErr &/*e*/) {
      //             std::cout << "Error @ " << x << " " << y << "\n";
      //             std::cout << "Exception caught: " << e.what() << "\n";
      //             std::cout << "PRERHS: " << pre_rhs << "\n";
      //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
      //             std::cout << "RHS: " << rhs << "\n";
      //             std::cout << "LHS: " << lhs << "\n\n";
      //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
      //             exit(0);
    }
  
    d += lhs; // update parameter
    iter ++;
  }
  return d;
} 








