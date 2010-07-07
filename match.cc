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


void UpdateMatchingParams(vector<vector<LOLAShot> > trackPts, string DRGFilename, ModelParams modelParams,GlobalParams globalParams)
{

   DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
   
   int k = 0;

   vector<Vector3> imgPts;
   imgPts = GetTrackPtsFromImage(trackPts[k],DRGFilename);
   
   ImageView<float> x_deriv = derivative_filter(DRG, 1, 0);
   ImageView<float> y_deriv = derivative_filter(DRG, 0, 1);
   
   Vector<float,6> d;//defines the affine transform
   d(0) = 1.0;
   d(4) = 1.0;
   Matrix<float,6,6> rhs;
   Vector<float,6> lhs;

   int ii, jj, x_base, y_base;
   //x_base and y_base are the initial 2D positions of the image points
   x_base = imgPts[0][1];
   y_base = imgPts[0][2];
   float xx = x_base + d[0] * ii + d[1] * jj + d[2];
   float yy = y_base + d[3] * ii + d[4] * jj + d[5];
   
  
   InterpolationView<EdgeExtensionView<DiskImageView<PixelMask<PixelGray<uint8> > > , ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
                interpolate(DRG, BilinearInterpolation(), ZeroEdgeExtension());

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
   d += lhs;
}










