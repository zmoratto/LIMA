// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace std;


//computes the translation between the foreground and background pixels
void ComputeDEMTranslation(vector<Vector3> featureArray, vector<Vector3> matchArray, Vector3 &translation );

void PrintMatrix(Matrix<float, 3, 3> &A);

//compute the matching error vector and overall average error 
float ComputeMatchingError(vector<Vector3> featureArray, vector<Vector3> matchArray, vector<float> &errorArray);

void ComputeDEMRotation(vector<Vector3> featureArray, vector<Vector3> matchArray, 
                        Vector3 translation, Matrix<float, 3, 3> &rotation);

//applies a 3D rotation and transform to a DEM
void  TransformFeatures(vector<Vector3> &featureArray, Vector3 translation, Matrix<float,3,3> rotation);



//compute a set of features from the fore image.
//TODO: copy featureVec in cartesian coordinates
template <class ViewT>
vector<Vector3> 
GetFeatures(ImageViewBase<ViewT> const& foreImg, GeoReference const &foreGeo,
            ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo)
{
 
 int verStep = 8;
 int horStep = 8;
 vector<Vector3> featureArray;
 
 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
  
 for (int j = 0; j < foreImg.impl().rows(); j=j+verStep){
     for (int i = 0; i < foreImg.impl().cols(); i=i+horStep){
       //determine the location of the corresponding background point

       if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != 0)  { 

	 Vector2 forePix(i,j);
	 Vector2 fore_lonlat = foreGeo.pixel_to_lonlat(forePix);
	 fore_lonlat(0) = fore_lonlat(0)-180;
       
         Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));
         //Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat3);
      
         featureArray.push_back(fore_lonlat3);
       }
     }
  }

 return featureArray;

};
//compute a set of features from the fore image.
//TODO: copy featureVec in cartesian coordinates
template <class ViewT>
vector<Vector3> 
GetFeaturesXYZ(ImageViewBase<ViewT> const& foreImg, GeoReference const &foreGeo,
               ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo)
{
 
 int verStep = 8;
 int horStep = 8;
 vector<Vector3> featureArray;
 
 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
  
 for (int j = 0; j < foreImg.impl().rows(); j=j+verStep){
     for (int i = 0; i < foreImg.impl().cols(); i=i+horStep){
       //determine the location of the corresponding background point

       if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != 0)  { 

	 Vector2 forePix(i,j);
	 Vector2 fore_lonlat = foreGeo.pixel_to_lonlat(forePix);
	 fore_lonlat(0) = fore_lonlat(0)-180;
       
         Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));
         Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat3);
      
         featureArray.push_back(fore_xyz);
       }
     }
  }

 return featureArray;

};
//find the closest background points to the fore points.
//TODO: copy matchingArray in cartesian coordinates
template <class ViewT>
void
FindMatches(vector<Vector3> featureArray, ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo, vector<Vector3>& matchArray )
{
 
 int matchWindowHalfWidth = 5;
 int matchWindowHalfHeight = 5;
 float usgs_2_lonlat = 180/(3.14159265*3396190);

 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 
 for (int i = 0; i < featureArray.size(); i++){
      
     Vector2 fore_lonlat;
     fore_lonlat(0) = featureArray[i](0); 
     fore_lonlat(1) = featureArray[i](1);
     Vector2 back_lonlat = fore_lonlat/usgs_2_lonlat;
     Vector2 backPix = backGeo.lonlat_to_pixel(back_lonlat);	 
     float x = backPix(0);
     float y = backPix(1);
     Vector3 back_lonlat3;
     
     //TO DO:revert fore from cartesian to spherical coordinates 
     //Vector3 back_lonlat3(fore_lonlat(0), fore_lonlat(1), (interpBackImg.impl())(x,y));

     //this can be commented out - START
     float minDistance = 10000000000.0;
     for (int k =  y-matchWindowHalfHeight; k < y + matchWindowHalfHeight; k++){
        for (int l = x - matchWindowHalfWidth; l < x + matchWindowHalfWidth; l++){
 
        
          Vector2 backPix(l,k);
	  Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
          Vector3 t_back_lonlat3;
          
          t_back_lonlat3(0) = back_lonlat(0)*usgs_2_lonlat; 
          t_back_lonlat3(1) = back_lonlat(1)*usgs_2_lonlat; 
          t_back_lonlat3(2) = interpBackImg(l,k);
           
          float distance1 = t_back_lonlat3(0) - featureArray[i](0);
          float distance2 = t_back_lonlat3(1) - featureArray[i](1);
          float distance3 = t_back_lonlat3(2) - featureArray[i](2);
          float distance = sqrt(distance1*distance1 + distance2*distance2 + distance3*distance3);   
          if (distance < minDistance){
	    minDistance = distance;
            back_lonlat3 = t_back_lonlat3;
          }
  
	}
     }
     //this can be commented out - END
     
     matchArray[i]=back_lonlat3;
 }
 
};


//find the closest background points to the fore points.
//TODO: copy matchingArray in cartesian coordinates
template <class ViewT>
void
FindMatchesXYZ(vector<Vector3> featureArray, ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo, 
               GeoReference const &foreGeo, vector<Vector3>& matchArray )
{
 
 int matchWindowHalfWidth = 5;
 int matchWindowHalfHeight = 5;
 float usgs_2_lonlat = 180/(3.14159265*3396190);

 ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 
 for (int i = 0; i < featureArray.size(); i++){

     //TO DO:revert fore from cartesian to spherical coordinates
     Vector3 fore_lonlat3 = foreGeo.datum().cartesian_to_geodetic(featureArray[i]);
     Vector2 fore_lonlat;
     fore_lonlat(0) = fore_lonlat3(0); 
     fore_lonlat(1) = fore_lonlat3(1);
     
     Vector2 back_lonlat = fore_lonlat/usgs_2_lonlat;
     Vector2 backPix = backGeo.lonlat_to_pixel(back_lonlat);	 
     float x = backPix(0);
     float y = backPix(1);
     Vector3 back_lonlat3;
       
     //this can be commented out - START
     float minDistance = 10000000000.0;
     for (int k =  y-matchWindowHalfHeight; k < y + matchWindowHalfHeight; k++){
        for (int l = x - matchWindowHalfWidth; l < x + matchWindowHalfWidth; l++){
 
        
          Vector2 backPix(l,k);
	  Vector2 back_lonlat = backGeo.pixel_to_lonlat(backPix);
          Vector3 t_back_lonlat3;
          
          t_back_lonlat3(0) = back_lonlat(0)*usgs_2_lonlat; 
          t_back_lonlat3(1) = back_lonlat(1)*usgs_2_lonlat; 
          t_back_lonlat3(2) = interpBackImg(l,k);
           
          float distance1 = t_back_lonlat3(0) - featureArray[i](0);
          float distance2 = t_back_lonlat3(1) - featureArray[i](1);
          float distance3 = t_back_lonlat3(2) - featureArray[i](2);
          float distance = sqrt(distance1*distance1 + distance2*distance2 + distance3*distance3);   
          if (distance < minDistance){
	    minDistance = distance;
            back_lonlat3 = t_back_lonlat3;
          }
  
	}
     }
     //this can be commented out - END
     
     matchArray[i]=back_lonlat3;
 }
 
};

#if 0
//applies a 3D rotation and transform to a DEM
template <class ViewT>
void 
TransformDEM(ImageViewBase<ViewT> const& dem, GeoReference const &geo,
                  Vector3 translation, Matrix<float,3,3> rotation)
{
   ImageView<PixelGray<float> > transfDEM(dem.cols(), dem.rows());

   for (int j = 0; j < dem.impl().rows(); j++){
     for (int i = 0; i < dem.impl().cols(); i++){
       
       //determine the location of the corresponding background point

       if ((isnan(dem.impl()(i,j))!=FP_NAN) && (dem.impl()(i,j)) != 0)  { 

	 Vector2 initPix(i,j);
	 Vector2 init_lonlat = geo.pixel_to_lonlat(initPix);
         Vector3 init_lonlat3(init_lonlat(0), init_lonlat(1), dem.impl()(i,j));

         Matrix<float,3,3> R;//  = cvCreateMat(3,3,CV_32FC1);
         Vector<float,3>  i_lonlat3;//  = cvCreateMat(3,1,CV_32FC1);
         Vector<float,3> t_lonlat3;//  = cvCreateMat(3,1,CV_32FC1);
         Vector<float,3> T;//  = cvCreateMat(3,1,CV_32FC1);
	 /*         
         CvMat* R  = cvCreateMat(3,3,CV_32FC1);
         CvMat* i_lonlat3  = cvCreateMat(3,1,CV_32FC1);
         CvMat* t_lonlat3  = cvCreateMat(3,1,CV_32FC1);
         CvMat* T  = cvCreateMat(3,1,CV_32FC1);
         */
     
         for (int k = 0; k < 3; k++){
	   T[k] = translation(k);
           i_lonlat3[k]=init_lonlat3(k);
           for (int l = 0; l < 3; l++){
	     R[k][i] = rotation(k,i);
           }
	 }
	 
         /*
         for (int k = 0; k < 3; k++){
	   T->data.fl[k] = translation(k);
           i_lonlat3->data.fl[k]=init_lonlat3(k);
           for (int l = 0; l < 3; l++){
	     R->data.fl[k*3+i] = rotation(k,i);
           }
	 }
	 */

	 t_lonlat3=R*i_lonlat3+T;
         //cvMatMulAdd(R, i_lonlat3, T, t_lonlat3);

         Vector2 transf_lonlat(t_lonlat3[0], t_lonlat3[1]);
         Vector2 transPix = geo.lonlat_to_pixel(transf_lonlat);
         transfDEM(transPix(0), transPix(1)) = t_lonlat3[2];//new transformed value
	 /*
         Vector2 transf_lonlat(t_lonlat3->data.fl[0], t_lonlat3->data.fl[1]);
         Vector2 transPix = geo.lonlat_to_pixel(transf_lonlat);
         transfDEM(transPix(0), transPix(1)) = t_lonlat3->data.fl[2];//new transformed value;
	 */
       }
     }
   }
}
# endif
/*
//computes the translation between the foreground and background pixels
template <class ViewT>
Vector3
ComputeDEMTranslation(ImageViewBase<ViewT> const& foreImg, GeoReference const &foreGeo,
                      ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo)
{

  float usgs_2_lonlat = 180/(3.14159265*3396190);

 
  ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 
  Vector3 translation;
  translation[0] = 0;
  translation[1] = 0;
  translation[2] = 0;

  int count = 0;
  for (int j = 0; j < foreImg.impl().rows(); j++){
     for (int i = 0; i < foreImg.impl().cols(); i++){
       //determine the location of the corresponding background point

       if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != 0)  { 

	 Vector2 forePix(i,j);
	 Vector2 fore_lonlat = foreGeo.pixel_to_lonlat(forePix);
	 fore_lonlat(0) = fore_lonlat(0)-180;
       
	 //Vector2 backPix = backGeo.lonlat_to_pixel(fore_lonlat);
        

	 Vector2 back_lonlat = fore_lonlat/usgs_2_lonlat;

       
	 //if (back_lonlat(0)< 0){ 
	 //    back_lonlat(0) = 180+back_lonlat(0);
	 //}
	 Vector2 backPix = backGeo.lonlat_to_pixel(back_lonlat);
       
	 float x = backPix(0);
	 float y = backPix(1);
         
         //Vector3 back_longlat3(back_lonlat(0), back_lonlat(1), (interpBackImg.impl())(x,y));
         //Vector3 back_xyz = backGeo.datum().geodetic_to_cartesian(back_longlat3);
         Vector3 back_lonlat3(fore_lonlat(0), fore_lonlat(1), (interpBackImg.impl())(x,y));
         Vector3 back_xyz = foreGeo.datum().geodetic_to_cartesian(back_lonlat3);

         Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));
         Vector3 fore_xyz = foreGeo.datum().geodetic_to_cartesian(fore_lonlat3);

	 //compute the local translation
         translation[0] = translation[0] + fore_lonlat3(0)-back_lonlat3(0);
	 translation[1] = translation[1] + fore_lonlat3(1)-back_lonlat3(1);
         translation[2] = translation[2] + fore_lonlat3(2)-back_lonlat3(2);
       
         count++;
       }
     }
  }
  translation[0] = translation[0]/count;
  translation[1] = translation[1]/count;
  translation[2] = translation[2]/count;
  return translation;

}
*/
/*
//computes the translation between the foreground and background pixels
template <class ViewT>
Matrix<float,3,3> 
ComputeDEMRotation(ImageViewBase<ViewT> const& foreImg, GeoReference const &foreGeo,
                   ImageViewBase<ViewT> const& backImg, GeoReference const &backGeo, 
                   Vector3 translation)
{

  float usgs_2_lonlat = 180/(3.14159265*3396190);

 
  ImageViewRef<typename ViewT::pixel_type>  interpBackImg = interpolate(edge_extend(backImg.impl(),
                                                                           ConstantEdgeExtension()),
                                                                           BilinearInterpolation());
 

  Matrix<float, 3, 3> rotation;
  Matrix<float,3,3> A;
  Matrix<float,3,3> U;
  Vector<float,3>   s;
  Matrix<float,3,3> V;

  
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      A[0][0] = 0.0;
    }
  }

  //#if 0
  int count = 0;
  //for each fore image point 
  for (int j = 0; j < foreImg.impl().rows(); j++){
     for (int i = 0; i < foreImg.impl().cols(); i++){

       //determine the location of the corresponding background point
       if ((isnan(foreImg.impl()(i,j)[0])!=FP_NAN) && (foreImg.impl()(i,j)[0]) != 0)  { 
	 Vector2 forePix(i,j);
	 Vector2 fore_lonlat = foreGeo.pixel_to_lonlat(forePix);
         //Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));

	 //Vector2 backPix = backGeo.lonlat_to_pixel(fore_lonlat);

	 fore_lonlat(0) = fore_lonlat(0)-180;
	 Vector2 back_lonlat = fore_lonlat/usgs_2_lonlat;
	 //if (back_lonlat(0)< 0){ 
	 //    back_lonlat(0) = 180+back_lonlat(0);
	 //}
	 Vector2 backPix = backGeo.lonlat_to_pixel(back_lonlat);
	 float x = backPix(0);
	 float y = backPix(1);
         Vector3 back_lonlat3(fore_lonlat(0), fore_lonlat(1), (interpBackImg.impl())(x,y));
         
         Vector3 fore_lonlat3(fore_lonlat(0), fore_lonlat(1), (foreImg.impl())(i,j));
	
         //TODO: update the translation vector as well
       
         A[0][0] = A[0][0] + back_lonlat3[0]*(fore_lonlat3[0] - translation[0]);
         A[0][1] = A[0][1] + back_lonlat3[1]*(fore_lonlat3[0] - translation[0]);
	 A[0][2] = A[0][2] + back_lonlat3[2]*(fore_lonlat3[0] - translation[0]);
	 A[1][0] = A[1][0] + back_lonlat3[0]*(fore_lonlat3[1] - translation[1]);
	 A[1][1] = A[1][1] + back_lonlat3[1]*(fore_lonlat3[1] - translation[1]);
	 A[1][2] = A[1][2] + back_lonlat3[2]*(fore_lonlat3[1] - translation[1]);
	 A[2][0] = A[2][0] + back_lonlat3[0]*(fore_lonlat3[2] - translation[2]);
	 A[2][1] = A[2][1] + back_lonlat3[1]*(fore_lonlat3[2] - translation[2]);
	 A[2][2] = A[2][2] + back_lonlat3[2]*(fore_lonlat3[2] - translation[2]);
	 
       }
     }
  }
  

  printf("correlation: A\n");
  PrintMatrix(A);
  
  //extract the rotation matrix
  //void svd( AMatrixT const& A, UMatrixT &U, SingularValuesT &s, VTMatrixT &VT );
  svd(A, U, s, V);
  

  printf("A = UWV^T\n");
  printf("W\n");
  PrintMatrix(A);
 
  printf("U\n");
  PrintMatrix(U);
 
  printf("V\n");
  PrintMatrix(V);

  printf("SVD done\n");

  Matrix<float,3,3> VT = transpose(V);
  PrintMatrix(VT);

  rotation = U*VT;

  printf("R\n");
  PrintMatrix(rotation);

  return rotation;

}
*/


