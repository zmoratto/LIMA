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
#include "tracks.h"
#include "weights.h"

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


void UpdateMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
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



  cout << "UMP: ComputeScaleFactor..." << endl;
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);


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






  int i_access, j_access;
  int ii, jj;
  int i_Center, j_Center;
  vector<int> center_ij;
  float x_base, y_base;
  float xx,yy;
  int row_max, col_max;

  //working with the entire CVS file significant errors where seen at index 126 - set index to this by hand to debugg this test case
 /*
  1. we see huge problems from the get go - massive lhs numbers, where does these come from in comparison
  2. need to explore & instrument the code further
  3. we have some stable/unstable pairs to explore (106 stable, 109 not) & ( 112 stable, 115 not)
 
  A couple of thoughts - we rezero lhs & rhs every iteration this should not be the problem.  Could the problem simple be a scaling issue where the size of the error terms on the right & left hand side - caused by the great number of points - causes numerical stability issues for solving the matrix?

  Try creating a seperate rhs_1k where you divide each element by 1000 directly before you attemt the inversion - what happens?  No relief, not surprising.

*/

  //the loop will start here
  //for (int index = 0; index < d2_array.size(); index++){
  for(int index = 112; index < 117; index +=3){

    //the following is for affine (NOT perspective) transforms
    Vector<float,6> d;//defines the affine transform
    d = d2_array[index];
    printf("d = [ %f, %f, %f, %f %f %f]\n",d[0],d[1],d[2],d[3],d[4],d[5]);

    const int tf_size = 6;//transform size
    Matrix<float,tf_size,tf_size> rhs;
    Matrix<float,tf_size,tf_size> rhs_1k;
    Vector<float,tf_size> lhs;
    Vector<float,tf_size> lhs_1k;
    print_rhs(rhs);

    //get row_max, col_max - improve this!
    cout << "UMP: max pixel locations..." << endl; 
    row_max = x_deriv.rows();
    col_max = x_deriv.cols();
    printf("row_max = %d, col_max = %d\n",row_max,col_max);

    int iter=0;
    cout << "UMP: index = "<< index << endl;
    cout << "UMP: grad descend loop ..." << endl; 

    // Calculate center points of the images
    cout << "UMP: access size, rows: " << row_max << ", cols:"<< col_max << endl;
    int i_C =  row_max/2;
    int j_C = col_max/2;
    cout << "Center at: (" << i_C << ", " << j_C << ")" << endl;

    int iA = 0;
    int jA = 0;

    // std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif"; 
    std::string filenameNoPath = sufix_from_filename(DRGFilename);
    std::string save_name_file = "../results" + prefix_less3_from_filename(filenameNoPath) + "lima_results"  +".txt";  
    std::string d_final_filename = "../results" + prefix_less3_from_filename(filenameNoPath) + "d_final" + ".txt";  

    FILE* sFile;
    sFile = fopen(save_name_file.c_str(),"w"); 

    FILE *d_FILE;
    d_FILE = fopen(d_final_filename.c_str(),"w"); // write the final result to the d_File at the end of this program

    float g_error = 0.0;
    float w_ght = 1.0; // folks you know what this'll be...! (i.e. weighted least squares) 
    

    //print out the max 
    while( iter <= numMaxIter){ //gradient descent => optimal transform

      int num_valid = 0;

      for (int k = 0; k < trackPts.size(); k++){

        for (int i = 0; i < trackPts[k].size(); i++){

          if ((trackPts[k][i].valid ==1) && (trackPts[k][i].reflectance !=0)){
           
            w_ght = trackPts[k][i].weight_prd;
           // if( (k==0) & (i==0)  ){
           //   printf("\n\n trackPts[%d][%d].weight_prd = %f =? w_ght = %f \n\n", k, i, trackPts[k][i].weight_prd, w_ght);
           // }
            w_ght = 1.0;
           

            num_valid += 1;
            // lola pixel coordinates

            for (int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){
              if (trackPts[k][i].LOLAPt[j].s == 3){//center point of a valid shot
                x_base = trackPts[k][i].imgPt[j].x;
                y_base = trackPts[k][i].imgPt[j].y;
              }   
            }
            //printf("x_base = %f, y_base = %f\n", x_base, y_base);
            // image coordinates under the current transform
            iA = (int) floor(d[0]*x_base + d[1]*y_base + d[2]);
            jA = (int) floor(d[3]*x_base + d[4]*y_base + d[5]);

            //printf("iA = %d, jA = %d\n", iA, jA);

            // caluculate ii & jj relative to the image center
            ii = iA - i_C;
            jj = jA - j_C;

            //let's chanck( ii, jj) are changing with the       
            if (num_valid < 0){
              // will the point be accepted into this round of the computation?
              int accept_comp = 0;
              if( ( iA >= 0) && ( iA <= row_max) && ( jA >= 0) && ( jA <= col_max)){
                accept_comp = 1; 
              }
              printf("iter = %d, num_valid = %d, pixel( %d) = ( %d, %d), accept_comp = %d\n",iter,num_valid,i,iA,jA,accept_comp);
            }
            // check (ii,jj) are inside the image!
            if( ( iA >= 0) && ( iA <= row_max) && ( jA >= 0) && ( jA <= col_max)){

              //initialize constants
              float I_x_sqr, I_x_I_y, I_y_sqr; 

              float I_e_val = interpDRG(jA,iA) - scaleFactor*trackPts[k][i].reflectance;
              //printf("I_e_val = %f\n", I_e_val);
              g_error += abs(I_e_val);
              float I_y_val, I_x_val;

              //calculate numerical dirivatives (ii,jj).
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
        }// end of for loop over i
      }//end of loop over k

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

      //copy over rhs_1k to see the outcome
      for(int i_rhs=0; i_rhs < tf_size; i_rhs++ ){
        for(int j_rhs=0;j_rhs < tf_size; j_rhs ++ ){
          rhs_1k( i_rhs, j_rhs) = rhs( i_rhs, j_rhs)/1000;
        }
        //lhs_1k(i_rhs) = lhs(i_rhs)/1000;
      }

      for( int i_lhs1K = 0; i_lhs1K < lhs_1k.size() ; i_lhs1K++){
        lhs_1k(i_lhs1K) = lhs(i_lhs1K)/1000;
      }


        try {
        solve_symmetric_nocopy(rhs,lhs);
        solve_symmetric_nocopy(rhs_1k,lhs_1k);
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


      printf("Checking if we can pull a linear factor of 1000 of all lhs & rhs enteries to achieve numerical stability - I sure hope we can!\n");
        //Make your changes here then print out and view the results! 

      //write out 
      fprintf(sFile,"iter= %d g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n",iter, g_error, d(0), d(1), d(2), d(3), d(4), d(5));
      //fclose(sFile);

      //see if you can print out this final problem
      printf("\n\n\nAre we hitting numeric difficulties?\n");
      printf("Numeric Work: iter %d: lhs_1k = [ %f, %f, %f, %f, %f, %f]\n",iter, lhs_1k[0], lhs_1k[1], lhs_1k[2], lhs_1k[3], lhs_1k[4], lhs_1k[5]);
      printf("Numeric Work: printing the rhs_1K to see what happens\n ");
      print_rhs(rhs_1k);
      printf("Above is numeric output\n\n\n ");

      error_array[index] = g_error;
      final_d_array[index] = d;

      g_error = 0.0;    

      cout << "UMP: Just written to latest... "<< endl;

      d += lhs; // update parameter - should this be d = lhs?

      printf("B. iter %d: d = [ %f, %f, %f, %f %f %f]\n",iter,d[0],d[1],d[2],d[3],d[4],d[5]);
      iter ++;

      //reset rhs & lhs before next iter
      for(int i_RHS = 0; i_RHS < tf_size; i_RHS++){
        for(int j_RHS = 0; j_RHS < tf_size; j_RHS++){
          rhs(i_RHS,j_RHS) = 0.0;
        }
        lhs(i_RHS) = 0.0;
      }
    }
    
    fclose(sFile);
    // here - write final 'd' 
    d_FILE = fopen(d_final_filename.c_str(),"w"); // write the final result to the d_File at the end of this program
    fprintf(d_FILE,"d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", d(0), d(1), d(2), d(3), d(4), d(5));
    fclose(d_FILE);

  }//index loop ends here

} 

//OLD FUNCTIONS BELOW THIS LINE
//============================================================================
#if 0
void UpdateMatchingParamsOld(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
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

  vector<Vector3> imgPts;

  int num_allImgPts = 0;
  for(int k = 0; k<trackPts.size();k++){
    num_allImgPts += trackPts[k].size();
  }
  imgPts.resize(num_allImgPts);  

  int index= 0;
  for(int k = 0; k < trackPts.size();k++){
    for(int i = 0; i < trackPts[k].size(); i++){ 
      if (trackPts[k][i].valid != 1){
        imgPts[index][0] = -1;
        imgPts[index][1] = -1;
        imgPts[index][2] = -1;
      }
      else{
        for (int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){

          if (trackPts[k][i].LOLAPt[j].s == 3){//center point of a valid shot

            imgPts[index][0] = trackPts[k][i].imgPt[j].val;
            imgPts[index][1] = trackPts[k][i].imgPt[j].x ;
            imgPts[index][2] = trackPts[k][i].imgPt[j].y;

          }   
        }
      }
      index++;
    }
  }

  //compute the synthetic image values
  //cout << "UMP: ComputeTrackReflectance..." << endl;
  //ComputeAllReflectance(trackPts, modelParams, globalParams);

  cout << "UMP: ComputeScaleFactor..." << endl;
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);

  cout << "UMP: ComputeSynthImgPts..." << endl;
  vector<float> synthImg = ComputeSyntImgPts(scaleFactor, trackPts);

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

  int i_access, j_access;
  int ii, jj;
  int i_Center, j_Center;
  vector<int> center_ij;
  float x_base, y_base;
  float xx,yy;
  int row_max, col_max;


  //the loop will start here
  for (int index = 0; index < d2_array.size(); index++){

    //the following is for affine (NOT perspective) transforms
    Vector<float,6> d;//defines the affine transform
    d = d2_array[index];
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
    cout << "UMP: index = "<< index << endl;
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

    std::string filenameNoPath = sufix_from_filename(DRGFilename);
    std::string save_name_file = "../results" + prefix_less3_from_filename(filenameNoPath) + "lima_results"  +".txt";  
    std::string d_final_filename = "../results" + prefix_less3_from_filename(filenameNoPath) + "d_final" + ".txt";  

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

            float I_e_val = interpDRG(jA,iA) - synthImg[i];
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
      fclose(sFile);


      error_array[index] = g_error;
      final_d_array[index] = d;

      g_error = 0.0;    

      cout << "UMP: Just written to latest... "<< endl;

      d += lhs; // update parameter - should this be d = lhs?
      printf("B. iter %d: d = [ %f, %f, %f, %f %f %f]\n",iter,d[0],d[1],d[2],d[3],d[4],d[5]);
      iter ++;

      //reset rhs & lhs before next iter
      for(int i_RHS = 0; i_RHS < tf_size; i_RHS++){
        for(int j_RHS = 0; j_RHS < tf_size; j_RHS++){
          rhs(i_RHS,j_RHS) = 0.0;
        }
        lhs(i_RHS) = 0.0;
      }
    }
    // here - write final 'd'
    fprintf(d_FILE,"d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", d(0), d(1), d(2), d(3), d(4), d(5));
    fclose(d_FILE);

  }//index loop ends here

} 

#endif


/*Notes on current behavior
  We are seeing two things we don't like:

  1. blatantly wrong answer on the first stage of the iteration

  2. (Solved) later iterations only appear add the output of the first stage to the previous result.  Current results point towards rhs not changing at all after iter 1!
  => does rhs change after iteration 1?

  No change in rhs after iter-1, we need to reset!

  3. Should we be using d = lhs or d += lhs as we are currently doing.  I need to think more about this. 

 */


