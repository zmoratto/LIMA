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

float ComputeScaleFactor(vector<float> allImgPts, vector<float> reflectance)
{
    float nominator =0.0;
    int numValidPts = 0;
    float scaleFactor = 1;

    for(int m = 0; m < reflectance.size(); m ++){
      if ((reflectance[m]!= -1) && (reflectance[m] !=0.0)){
          nominator += allImgPts[m]/reflectance[m];
          numValidPts += 1;
      }
    }

    if (numValidPts != 0){ 
       scaleFactor = nominator/numValidPts;
    }

    return scaleFactor;
}

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


Vector<float,6> UpdateMatchingParams(vector<vector<LOLAShot> > trackPts, string DRGFilename, ModelParams modelParams,GlobalParams globalParams)
{
   DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
   
   int k = 0;
   
   //get the true image points
   //cout << "UMP calling: GetTrackPtsFromImage" << endl;
   //cout << "rather sure this function doesn't exist here...!"<< endl;
   vector<Vector3> imgPts;
   imgPts = GetTrackPtsFromImage(trackPts[k],DRGFilename);
   
   //cout << "UMP finished: GetTrackPtsFRomImage " << endl; 

   //compute the synthetic image values
   cout << "UMP: ComputeTrackReflectance..." << endl;
  vector<float> reflectance = ComputeTrackReflectance(trackPts[k], modelParams, globalParams);
    cout << "UMP: ComputeScaleFactor..." << endl;
float scaleFactor = ComputeScaleFactor(imgPts, reflectance);
      //isn't this a bit self-serving?  The scale factor should be computed a priori or conditioned from a learned distribution.  Here, we are confusing our training & testing data.
    cout << "UMP: ComputeSynthImgPts..." << endl;
  vector<float> synthImg = ComputeSyntImgPts(scaleFactor, reflectance);
       
   cout << "UMP: calculate min & max pixel locations..."<< endl;    
   
   int i_min = 0;
   int i_max = DRG.rows();
   int j_min = 0;
   int j_max = DRG.cols();
   

   printf("x range from: %f -> %f \ny range from: %f -> %f \n",i_min,i_max,j_min,j_max);


   cout << "UMP: derivative_filter "<< endl; 
   ImageView<float> x_deriv = derivative_filter(DRG, 1, 0);
   ImageView<float> y_deriv = derivative_filter(DRG, 0, 1);

  //the following is for affine (NOT perspective) transforms

   Vector<float,6> d;//defines the affine transform
   d(0) = 1.0;
   d(4) = 1.0;
   Matrix<float,6,6> rhs;
   Vector<float,6> lhs;

   int ii, jj;
   float x_base, y_base;
   float xx,yy;
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
  while(iter < 20) // gradient descent loop to iterate to optimal transform
  {

   for (int i = 0; i < imgPts.size(); i++)
   {
  
        //
        if( (synthImg[i]!= -1) && (synthImg[i]!=0) )
        {
          // the above if statement will eventually encompas everything!
        
        // lola pixel coordinates
        x_base = imgPts[i][1];
        y_base = imgPts[i][2];
        
        // find image coordinates of the lola point based on the current transform
        ii = (int) floor(d[0]*x_base + d[1]*y_base + d[2]);
        
        jj = (int) floor(d[3]*x_base + d[4]*y_base + d[5]);
        

        // load up (x_base, y_base)
        // if(){} // check that the points is 1. valid & 2. !0

        float I_x_sqr, I_x_I_y, I_y_sqr; 
        float I_e_val = imgPts[i][0]-synthImg[i];
        float I_y_val, I_x_val;
        
        //calculate numerical dirivatives (ii,jj)... 
        printf("i = %d: (%d,%d)",ii,jj);
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
   
        }// end of if statement
   }// end of for loop over all data  points
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
   d += lhs; // update parameters
   
   // calculate stopping condition - currently at just 
    iter ++;
  }// end of while loop  updating parameters
 return d;
} //end of function












