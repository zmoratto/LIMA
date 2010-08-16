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

#include <stdio.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
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

template <class ViewT>
vector<Vector3>
GetAllPtsFromImage( vector<vector<LOLAShot > > trackPts,  ImageViewBase<ViewT> const& DRG, GeoReference const &DRGGeo)
{
 
  //typedef typename PixelChannelType<typename ViewT::pixel_type>::type channel_type;

  vector<Vector3> allImgPts;

  int num_allImgPts = 0;
  for(int k = 0; k<trackPts.size();k++)
  {
    num_allImgPts += trackPts[k].size();
  }

  allImgPts.resize(num_allImgPts);  
  vector<pointCloud> ptHere;

  //ImageViewRef<PixelMask<PixelGray<uint8> > >  interpDRG = interpolate(edge_extend(DRG.impl(),
  //      ConstantEdgeExtension()),
  //    BilinearInterpolation());
  /*
  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
                                                            ConstantEdgeExtension()),
                                                            BilinearInterpolation())
  */

  ImageViewRef<float> interpDRG;
  if ( IsMasked<typename ViewT::pixel_type>::value == 0 ) {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(DRG.impl(),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );

    cout << "NOT masked" <<endl;
  } else {
    interpDRG = pixel_cast<float>(interpolate(edge_extend(apply_mask(DRG.impl()),
							  ConstantEdgeExtension()),
					      BilinearInterpolation()) );
     cout << "MASKED" <<endl;
  }

  int i_tot = 0;
  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      ptHere = trackPts[k][i].LOLAPt;
      pointCloud centerPt  = GetPointFromIndex( ptHere, 3);
      pointCloud topPt     = GetPointFromIndex( ptHere, 2);
      pointCloud leftPt    = GetPointFromIndex( ptHere, 1);

      if((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1)){

        float lon = centerPt.coords[0];
        float lat = centerPt.coords[1];
        float rad = centerPt.coords[2];
     
        Vector2 DEM_lonlat(lon, lat);
        Vector2 DRG_pix = DRGGeo.lonlat_to_pixel(DEM_lonlat);

        //int x = (int)DRG_pix[0];
        //int y = (int)DRG_pix[1];
        //cout << "lon="<< lon << "lat=" << lat << endl; 
        
        float x = DRG_pix[0];
        float y = DRG_pix[1];
        //cout << "x="<< x << "y=" << y << endl; 
        //PixelMask<PixelGray<uint8> > DRGVal = interpDRG(x, y);
        //PixelGray<uint8> DRGVal = interpDRG(x, y);
        //Pixel<PixelGray<uint8> > DRGVal = interpDRG(x,y);


        //insert data
        //allImgPts[i_tot + i][0] = (float) DRGVal;
        //cout<< "interp = " << interpDRG(x,y) << endl;
        allImgPts[i_tot + i][0] = interpDRG(x, y);
        allImgPts[i_tot + i][1] = DRG_pix[0];
        allImgPts[i_tot + i][2] = DRG_pix[1];

        /* Debugg statements to check data is correctly copied
        //The incorrect cast from int -> float is root cause
        if(i%15==0)
        {
        printf("LOLA # %d @(lon = %f , lat = %f) = pixel(%d,%d)\nwith incorrect cast at pixel(%f,%f)\n",i,lon,lat,x,y,x,y);
        }
        //
        //
        if(i%30==0)
        {
        cout <<"Is DRG_pix fractional?" << endl;        
        printf("i = %d, DRG_pix = [%f,%f], pixel = (%d,%d)\n",i,DRG_pix[0],DRG_pix[1],x,y);
        }
         */
      }else {
        //write -1 to designate and invalid point
        allImgPts[i_tot + i][0] = -1;
        allImgPts[i_tot + i][1] = -1;
        allImgPts[i_tot + i][2] = -1;
      }
    }
    i_tot += trackPts[k].size();  
  }

  /* debugg statement: is data correctly copied?
     printf("GetAllPtsFromImg: print some pnts...\n");
     for(int i = 0; i< allImgPts.size(); i++)
     {
     if(i%150==0)
     {
     printf("allImgPts[%d][ %f, %f, %f]\n", i, allImgPts[i][0],  allImgPts[i][1], allImgPts[i][2]);
     }
     }
   */

  return allImgPts;
}


bool deriv_cached(string & input_name, string & cached_table){
  bool identified = false;
  string line;
  fstream myfile(cached_table.c_str());
  // does the name exist
  if(myfile.is_open()){
    while(! myfile.eof()){
      getline(myfile,line);
      cout << line << endl;
      if(line == input_name){
        identified = true;
        cout << "Fount it!" << endl; 
      }
    }
  }
  myfile.close();
  
  //if not add the name to the cached file
  if(!identified){
    // if we didn't find the file append the file to the document of saves
    ofstream file_change;
    file_change.open(cached_table.c_str(),ios::out | ios:: app);
    file_change << input_name << "\n";
    file_change.close();
  }
  

  return identified;
}



vector<float> ComputeAllReflectance( vector< vector<LOLAShot> >  allTracks, ModelParams modelParams, GlobalParams globalParams)
{
  vector<float> reflectance;
  vector<float> single_track_ref;
  int num_reflc = 0;

  for(int k = 0; k < allTracks.size(); k++ )
  {
    num_reflc += allTracks[k].size(); 
  }
  //cout << "num_reflc: " << num_reflc << endl;

  //reflectance.resize();
  //cout << "reflectance.size() should be: " << reflectance.size() << endl; 
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


void print_rhs( Matrix<float,6,6> rhs)
{
  printf("Printing out RHS:\n");

  printf("[ %f %f %f] [ %f %f %f]\n",rhs(0,0), rhs(0,1), rhs(0,2), rhs(0,3), rhs(0,4), rhs(0,5) );   
  printf("[ %f %f %f] [ %f %f %f]\n", rhs(1,0), rhs(1,1), rhs(1,2), rhs(1,3), rhs(1,4),rhs(1,5) ); 
  printf("[ %f %f %f] [ %f %f %f]\n\n", rhs(2,0), rhs(2,1),rhs(2,2), rhs(2,3), rhs(2,4), rhs(2,5) );

  printf("[ %f %f %f] [ %f %f %f]\n", rhs(3,0), rhs(3,1), rhs(3,2), rhs(3,3), rhs(3,4), rhs(3,5) );
  printf("[ %f %f %f] [ %f %f %f]\n", rhs(4,0), rhs(4,1), rhs(4,2), rhs(4,3), rhs(4,4), rhs(4,5) );
  printf("[ %f %f %f] [ %f %f %f]\n\n", rhs(5,0), rhs(5,1), rhs(5,2), rhs(5,3), rhs(5,4), rhs(5,5) );
}



void UpdateMatchingParams(vector<vector<LOLAShot> > trackPts, string DRGFilename,  
                     ModelParams modelParams, GlobalParams globalParams, int numMaxIter, 
                     vector<Vector<float, 6> >d2_array, vector<Vector<float, 6> >&final_d_array, 
                     vector<float> &error_array )
{
 
  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);
  

   ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
                                                                         ConstantEdgeExtension()),
  							    BilinearInterpolation());

  cout <<"Done interpolating the subsampled image"<<endl;
  

  //get the true image points
  cout << "UMP calling: GetAllPtsFromImage..." << endl;
  vector<Vector3> imgPts;

  imgPts = GetAllPtsFromImage(trackPts,  interpDRG, DRGGeo);

  //compute the synthetic image values
  cout << "UMP: ComputeTrackReflectance..." << endl;
  vector<float> reflectance = ComputeAllReflectance( trackPts, modelParams, globalParams);

  cout << "UMP: ComputeScaleFactor..." << endl;
  float scaleFactor = ComputeScaleFactor(imgPts, reflectance);
  //isn't this a bit self-serving?  The scale factor should be computed a priori or conditioned from a learned distribution.  Here, we confuse training & test data.

  cout << "UMP: ComputeSynthImgPts..." << endl;
  vector<float> synthImg = ComputeSyntImgPts(scaleFactor, reflectance);


  cout << "UMP: derivative_filter..."<< endl; 
  
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DRGFilename);

  std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_y_deriv.tif";

  cout << "sufix_from_filename(DRGFilename): "<< sufix_from_filename(DRGFilename) << endl;
  cout << "UPM: temp = " << temp << endl;
  cout << "prefix_less3_from_filename(temp): "<< prefix_less3_from_filename(temp) << endl;

  cout<<xDerivFilename<<endl;

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "UMP: trying to compute x_derivative" << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
                                temp.format(), Vector2i(512,512) );
    cout << "UMP: trying to write x_derivative" << endl;
    block_write_image(rsrc, temp,
                      TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
  }
  cout << "Trying to read in x_deriv( xDerivFilename ) " << endl;
  
  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
                                temp.format(), Vector2i(512,512) );
    block_write_image(rsrc, temp,
                      TerminalProgressCallback("asp", "Derivative:") );
  }
  DiskImageView<float> y_deriv( yDerivFilename );

/*
  // single threaded version - should work better on Dave's machine
  // we have strings xDerivFilename & yDerivFilename
  cout << "UMP: STARTING NEW DERIVATIVE CODE!!" << endl;
  ImageView<float> x_deriv;
  ImageView<float> y_deriv;
  
  string hold_file = "derivate_cache_names.txt";

  if(deriv_cached( xDerivFilename, hold_file)){
    cout << "UMP: found "<< xDerivFilename << endl;
    cout << "UMP: loading " << xDerivFilename <<  endl;
       
    read_image( x_deriv,  xDerivFilename);

  }else{
    //append cache list - accomplished in deriv_cache   
    cout << "UMP: did not find " << xDerivFilename << endl;
    cout << "UMP: Compute & Cache " << xDerivFilename  << endl;
    
    x_deriv = derivative_filter(pixel_cast_rescale<float>( DRG), 1, 0);
    write_image( xDerivFilename, x_deriv, TerminalProgressCallback("vw", "saving x_deriv: "));
  }

  if(deriv_cached( yDerivFilename, hold_file)){
    cout << "UMP: found "<< yDerivFilename << endl;
    cout << "UMP: Compute & Cache " << xDerivFilename << endl;
    
    read_image( y_deriv, yDerivFilename);

  }else{ 
    //append cache list - accomplished in deriv_cached
    cout << "UMP: did not find " << yDerivFilename << endl;
    cout << "UMP: Compute & Cache " << yDerivFilename << endl;
  
    y_deriv = derivative_filter(pixel_cast_rescale<float>( DRG), 0, 1);
    write_image( yDerivFilename, y_deriv, TerminalProgressCallback("vw", "saving y_deriv: "));
    
  }
 
*/
/* // old & incomplete derivative caching code - used as a model for the single threaded version that appears above.
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
*/
  int i_access, j_access;
  int ii, jj;
  int i_Center, j_Center;
  vector<int> center_ij;
  float x_base, y_base;
  float xx,yy;
  int row_max, col_max;


  //the loop will start here
 
  //the following is for affine (NOT perspective) transforms
  Vector<float,6> d;//defines the affine transform
  d = d2_array[0];
  printf("d = [ %f, %f, %f, %f %f %f]\n",d[0],d[1],d[2],d[3],d[4],d[5]);
 
  const int tf_size = 6;//transform size
  Matrix<float,tf_size,tf_size> rhs;
  Vector<float,tf_size> lhs;
  print_rhs(rhs);

  //get row_max, col_max - improve this!
  cout << "UMP: max pixel locations..." << endl; 
  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
  printf("row_max = %d, col_max = %d\n",row_max,col_max);

  int iter=0;
  cout << "UMP: grad descend loop ..." << endl; 

  // Calculate center points of the images
  cout << "UMP: access size, rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C =  row_max/2;
  int j_C = col_max/2;
  cout << "Center at: (" << i_C << ", " << j_C << ")" << endl;

  int iA = 0;
  int jA = 0;
  
  //Create disk view resources -access image pnts by pixel
  
  //open stuff up...skip automatic naming!
  
  //string result_file = ;
 
  // std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  
  std::string temp_name_sNow = sufix_from_filename(DRGFilename);
  std::string save_name_file = "../results" + prefix_less3_from_filename(temp_name_sNow) + "lima_results.txt";  
  std::string d_final_filename = "../results" + prefix_less3_from_filename(temp_name_sNow) + "d_final.txt";  

  FILE* sFile;
  sFile = fopen(save_name_file.c_str(),"w"); 
  
  FILE *d_FILE;
  d_FILE = fopen(d_final_filename.c_str(),"w"); // write the final result to the d_File at the end of this program
 
  float g_error = 0.0; 
  while( iter <= numMaxIter) //gradient descent => optimal transform
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

        // calculate access 

        // image coordinates under the current transform
        iA = (int) floor(d[0]*x_base + d[1]*y_base + d[2]);
        jA = (int) floor(d[3]*x_base + d[4]*y_base + d[5]);

        /* pulled out of 'for' loop - constant reference frame for d += lhs
           ii = (int) floor(d[0]*x_base + d[1]*y_base + d[2]);
           jj = (int) floor(d[3]*x_base + d[4]*y_base + d[5]);
         */

        // caluculate ii & jj
        ii = iA - i_C;
        jj = jA - j_C;

        /* debugg statement: checks data transcription    
           if( num_valid % 50 == 0){
        //printf("LOLA pnt: %d of intensity %f, at (%f,%f) transformed to (%d,%d)\n",i,imgPts[i][0],x_base,y_base,ii,jj); 
        printf("LOLA pnt: %d of intensity %f, at (%f,%f) transformed to (%d,%d)\n",i,imgPts[i][0],imgPts[i][1],imgPts[i][2],ii,jj);
        }
         */ 

        //let's chanck( ii, jj) are changing with the       
        if(num_valid < 0)
        {
          // will the point be accepted into this round of the computation?
          int accept_comp = 0;
          if( ( iA >= 0) && ( iA <= row_max) && ( jA >= 0) && ( jA <= col_max)){
            accept_comp = 1; 
          }
        printf("test3\n");
          printf("inter = %d, num_valid = %d, pixel( %d) = ( %d, %d), accept_comp = %d\n",iter,num_valid,i,iA,jA,accept_comp);
        }
        // check (ii,jj) are inside the image!
        if( ( iA >= 0) && ( iA <= row_max) && ( jA >= 0) && ( jA <= col_max)){

          //initialize constants
          float I_x_sqr, I_x_I_y, I_y_sqr; 
          float I_e_val = imgPts[i][0]-synthImg[i];
          //above line is incorrect - must access pixel (jA,iA):
          //I_e_val = DRG_img_view_resourse(jA,iA) - synthImg[i];
          
          I_e_val = interpDRG(jA,iA) - synthImg[i];
          g_error += abs(I_e_val);
          float I_y_val, I_x_val;

          //calculate numerical dirivatives (ii,jj)... 
          I_x_val = x_deriv(jA,iA); 
          I_y_val = y_deriv(jA,iA); 

          //cout << I_x_val <<":" << I_y_val << endl;
          I_x_I_y = I_x_val*I_y_val;        
          I_x_sqr = I_x_val*I_x_val;
          I_y_sqr = I_y_val*I_y_val;

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


          if( num_valid < 0){
            // print out stuff, questions to ask:
            //1. Why the zeros in the upper left?
            //2. What about the Inf/NaN in the lr

            printf("num_valid %d: I_x_val = %f,  I_x_sqr = %f\n",num_valid,I_x_val,I_x_sqr);
          }

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

    //print out previous and updates
    //printf("A. iter %d: d = [ %f, %f, %f, %f %f %f]\n",iter,d[0],d[1],d[2],d[3],d[4],d[5]);
    printf("A. iter %d: lhs = [ %f, %f, %f, %f, %f, %f]\n",iter, lhs[0], lhs[1], lhs[2], lhs[3], lhs[4], lhs[5]);
    print_rhs(rhs);
    
    //write out 
    fprintf(sFile,"iter= %d g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n",iter, g_error, d(0), d(1), d(2), d(3), d(4), d(5));
    
    error_array[0] = g_error;
    final_d_array[0] = d;

    g_error = 0.0;    

    cout << "UMP: Just written to latest... "<< endl;
    
    d += lhs; // update parameter - should this be d = lhs?
    printf("B. iter %d: d = [ %f, %f, %f, %f %f %f]\n",iter,d[0],d[1],d[2],d[3],d[4],d[5]);
    iter ++;

    //reset rhs & lhs before next iter
    for(int i_RHS = 0; i_RHS < tf_size; i_RHS++)
    {
      for(int j_RHS = 0; j_RHS < tf_size; j_RHS++)
      {
        rhs(i_RHS,j_RHS) = 0.0;
      }
      lhs(i_RHS) = 0.0;
    }

  }
  // here - write final 'd'
  fprintf(d_FILE,"d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", d(0), d(1), d(2), d(3), d(4), d(5));
  

  //loop will end here

} 



/*Notes on current behavior
  We are seeing two things we don't like:

  1. blatantly wrong answer on the first stage of the iteration

  2. (Solved) later iterations only appear add the output of the first stage to the previous result.  Current results point towards rhs not changing at all after iter 1!
  => does rhs change after iteration 1?

  No change in rhs after iter-1, we need to reset!

  3. Should we be using d = lhs or d += lhs as we are currently doing.  I need to think more about this. 

 */


//UNUSED FUNCTIONS BELOW
//========================================================================================================================================================================
/*
vector<int> pixel_center_LOLA_pnts( vector<Vector3> imgPts, Vector<float,6> d)
{
  vector<int> center_ij;
  center_ij.resize(2);
  float i_c = 0.0;
  float j_c = 0.0;
  int num_S = imgPts.size();
  //get average pixel location after current transform
  for(int i = 0; i < num_S; i++)
  {
    i_c += d[0]*imgPts[i][1] + d[1]*imgPts[i][2] + d[2];  
    j_c += d[3]*imgPts[i][1] + d[4]*imgPts[i][2] + d[5];
  }

  center_ij[0] = (int) floor( i_c / (float) num_S);
  center_ij[1] = (int) floor( j_c / (float) num_S);

  return center_ij;

}
*/
#if 0
vector<Vector3>  GetAllPtsFromImage(vector<vector<LOLAShot > > trackPts,string DRGFilename)
{ 
  //vector<Vector3> GetTrackPtsFromImage(vector<LOLAShot> trackPts, string DRGFilename)
  

  //DiskImageView<PixelMask<PixelGray<uint8> > >  DRG(DRGFilename);
  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
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

  //ImageViewRef<PixelMask<PixelGray<uint8> > >  interpDRG = interpolate(edge_extend(DRG.impl(),
  //      ConstantEdgeExtension()),
  //    BilinearInterpolation());

  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  int i_tot = 0;
  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){
      ptHere = trackPts[k][i].LOLAPt;
      pointCloud centerPt  = GetPointFromIndex( ptHere, 3);
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

        //PixelMask<PixelGray<uint8> > DRGVal = interpDRG(x, y);
        PixelGray<uint8> DRGVal = interpDRG(x, y);
        //Pixel<PixelGray<uint8> > DRGVal = interpDRG(x,y);

        //insert data
        allImgPts[i_tot + i][0] = (float) DRGVal;
        allImgPts[i_tot + i][1] = DRG_pix[0];
        allImgPts[i_tot + i][2] = DRG_pix[1];

        /* Debugg statements to check data is correctly copied
        //The incorrect cast from int -> float is root cause
        if(i%15==0)
        {
        printf("LOLA # %d @(lon = %f , lat = %f) = pixel(%d,%d)\nwith incorrect cast at pixel(%f,%f)\n",i,lon,lat,x,y,x,y);
        }
        //
        //
        if(i%30==0)
        {
        cout <<"Is DRG_pix fractional?" << endl;        
        printf("i = %d, DRG_pix = [%f,%f], pixel = (%d,%d)\n",i,DRG_pix[0],DRG_pix[1],x,y);
        }
         */
      }else {
        //write -1 to designate and invalid point
        allImgPts[i_tot + i][0] = -1;
        allImgPts[i_tot + i][1] = -1;
        allImgPts[i_tot + i][2] = -1;
      }
    }
    i_tot += trackPts[k].size();  
  }

  /* debugg statement: is data correctly copied?
     printf("GetAllPtsFromImg: print some pnts...\n");
     for(int i = 0; i< allImgPts.size(); i++)
     {
     if(i%150==0)
     {
     printf("allImgPts[%d][ %f, %f, %f]\n", i, allImgPts[i][0],  allImgPts[i][1], allImgPts[i][2]);
     }
     }
   */

  return allImgPts;
}

#endif
