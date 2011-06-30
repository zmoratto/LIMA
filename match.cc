// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <omp.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Photometry.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include <stdio.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
using namespace vw::photometry;

using namespace std;
#include <math.h>
#include "coregister.h"
#include "tracks.h"
#include "weights.h"


#define CHUNKSIZE   1


void SaveReportFile(vector<vector<LOLAShot> > &trackPts, vector<Vector<float, 6> >finalTransfArray,  
                      vector<float> errorArray, string reportFilename)
{
    FILE *fp = fopen(reportFilename.c_str(),"w");
    int numElements = errorArray.size();
    int bestResult = 0;
    float smallestError = 1000000000.0;

    for (int i = 0; i < numElements; i++){
       fprintf(fp,"index=%d d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f e=%f\n", 
	               i,
                       finalTransfArray[i](0), finalTransfArray[i](1), 
                       finalTransfArray[i](2), finalTransfArray[i](3),
                       finalTransfArray[i](4), finalTransfArray[i](5),
	               errorArray[i]);
        if ( (errorArray[i] < smallestError) && (errorArray[i] > 0) ){
          smallestError = errorArray[i]; 
          bestResult = i;
        }    
      
    }
    //print the best transform index and smallest error
    fprintf(fp, "bestTransformIndex = %d, smallest error = %f\n", bestResult, errorArray[bestResult]);
    //print the total number of tracks
    fprintf(fp, "numTotalTracks = %d\n", (int)trackPts.size());
    
    //print the number of features per track
    for (unsigned int ti = 0; ti < trackPts.size(); ti++){
      fprintf(fp,"trackID: %d ", ti);
      fprintf(fp,"numShots: %d ", (int)trackPts[ti].size());
      int numFeatures = 0;
      for (unsigned int si = 0; si < trackPts[ti].size(); si++){
	if ((trackPts[ti][si].featurePtRefl == 1.0) && (trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0)) { 
	  numFeatures++;  
	}
      }
      fprintf(fp,"numFeatures: %d\n", numFeatures);
   }


   fclose(fp);
}




/*
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
*/
void printRHS( Matrix<float,6,6> rhs,  int index, int iter)
{
  cout<<"index= "<<index<<", iter= "<<iter<<endl;
  printf("--------------------------------------\n");
  printf("RHS = \n");
  printf("[ %f %f %f | %f %f %f]\n",rhs(0,0), rhs(0,1), rhs(0,2), rhs(0,3), rhs(0,4), rhs(0,5) );   
  printf("[ %f %f %f | %f %f %f]\n", rhs(1,0), rhs(1,1), rhs(1,2), rhs(1,3), rhs(1,4),rhs(1,5) ); 
  printf("[ %f %f %f | %f %f %f]\n", rhs(2,0), rhs(2,1),rhs(2,2), rhs(2,3), rhs(2,4), rhs(2,5) );
  printf("[ %f %f %f | %f %f %f]\n", rhs(3,0), rhs(3,1), rhs(3,2), rhs(3,3), rhs(3,4), rhs(3,5) );
  printf("[ %f %f %f | %f %f %f]\n", rhs(4,0), rhs(4,1), rhs(4,2), rhs(4,3), rhs(4,4), rhs(4,5) );
  printf("[ %f %f %f | %f %f %f]\n", rhs(5,0), rhs(5,1), rhs(5,2), rhs(5,3), rhs(5,4), rhs(5,5) );
  printf("--------------------------------------\n\n");
}

void printLHS_Error( Vector<float,6> lhs, float error, int index, int iter)
{
 
  cout<<"index= "<<index<<", iter= "<<iter<<", error="<<error<<endl;
  printf("--------------------------------------\n");
  printf("LHS = [ %f, %f, %f, %f, %f, %f]\n", lhs[0], lhs[1], lhs[2], lhs[3], lhs[4], lhs[5]);
  printf("--------------------------------------\n\n");
}
/*
//This functions needs to be changed to generate multiple starting points for each 
//element of the affine transform 
void GenerateInitTransforms( vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings)
{
    initTransfArray.resize(settings.maxNumStarts);
    for (int i = 0; i < settings.maxNumStarts; i++){
      initTransfArray[i][0] = 1.0;
      initTransfArray[i][1] = 0.0;
      initTransfArray[i][2] = (i-settings.maxNumStarts/2)*2;
      initTransfArray[i][3] = 0.0;
      initTransfArray[i][4] = 1.0;
      initTransfArray[i][5] = 0.0;//(i-maxNumStarts/2)*25;
    }  
}
*/
//This functions needs to be changed to generate multiple starting points for each 
//element of the affine transform 
void GenerateInitTransforms(vector<Vector<float, 6> > &initTransfArray, CoregistrationParams settings)
{

  vector<vector<float> > params;
  params.resize(6);
 
  params[0].resize(1);
  params[0][0] = 1;  

  params[1].resize(1);
  params[1][0] = 0;

  params[2].resize(19);
  for (unsigned int i = 0; i < params[2].size(); i++){
    params[2][i] = (float)(-9.0+i);
  }
 
  params[3].resize(1);
  params[3][0] = 0;  

  params[4].resize(1);
  params[4][0] = 1;

  params[5].resize(19);
  for (unsigned int i = 0; i < params[5].size(); i++){
    params[5][i] = (float)(-9.0+i);
  }
 

  int maxNumStarts = 1;
  for( int i = 0; i < 6; i++){
    maxNumStarts = maxNumStarts*params[i].size();
  } 
  initTransfArray.resize(maxNumStarts);
  
  int index = 0;
  for (unsigned int i0 = 0; i0 < params[0].size(); i0++){
    for (unsigned int i1 = 0; i1 < params[1].size(); i1++){ 
      for (unsigned int i2 = 0; i2 < params[2].size(); i2++){
        for (unsigned int i3 = 0; i3 < params[3].size(); i3++){
          for (unsigned int i4 = 0; i4 < params[4].size(); i4++){
            for (unsigned int i5 = 0; i5 < params[5].size(); i5++){
	      initTransfArray[index][0] = params[0][i0];
	      initTransfArray[index][1] = params[1][i1];
	      initTransfArray[index][2] = params[2][i2];
	      initTransfArray[index][3] = params[3][i3];
	      initTransfArray[index][4] = params[4][i4];
	      initTransfArray[index][5] = params[5][i5];
              index++;
	    }
	  }
	}
      }
    }
  }
  /*
  for (int i = 0; i < maxNumStarts; i++){
    printf("%d:", i);
    for (int j = 0; j < 6; j++){
      printf("%f ", initTransfArray[i][j]);
    }
    printf("\n");
  }
  */
}

void GetBestTransform(vector<Vector<float, 6> > &finalTransfArray,  vector<float> &finalMatchingErrorArray, 
                      Vector<float, 6> &optimalTransf, float &optimalError)
{
    //this will be made a special function to make sure it is flexible wrt to the name of the parameters.
    int bestResult = 0;
    float smallestError = std::numeric_limits<float>::max();
    for (unsigned int index = 0; index < finalTransfArray.size(); index++){
      printf("refined %d: g_error= %f d[0]= %f d[1]= %f d[2]= %f d[3]= %f d[4]= %f d[5]= %f\n", 
	     index, finalMatchingErrorArray[index], 
	     finalTransfArray[index](0), finalTransfArray[index](1), 
	     finalTransfArray[index](2), finalTransfArray[index](3),
	     finalTransfArray[index](4), finalTransfArray[index](5));
      if  ( (finalMatchingErrorArray[index] < smallestError) && (finalMatchingErrorArray[index] > 0) ){
	smallestError = finalMatchingErrorArray[index]; 
	bestResult = index;
      }    
    } 
    cout<<"bestResult= "<<bestResult<<endl;

    //copy the best transform to optimalTransfArray
    optimalTransf = finalTransfArray[bestResult];
    optimalError  = finalMatchingErrorArray[bestResult];
}

//determines the best finalTransfArray from all initTransfArrays
//it assumes that each image is transformed by one affine transform
//in the future we can investigate the use of one affine transform per track. 
void InitMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename,  
			       ModelParams modelParams, CoregistrationParams coregistrationParams,  
			       vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> > &finalTransfArray, 
			       vector<float> &matchingErrorArray)
{
    unsigned int ti, si, li;
    int jA, iA;
    float I_e_val;
    Vector2 minmax = ComputeMinMaxValuesFromCub(cubFilename);

    vector<float> errorArray;
    errorArray.resize(initTransfArray.size());
    vector<int> numValidPts;
    numValidPts.resize(initTransfArray.size());    

    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
    double nodata_value = rsrc->nodata_read();
    cout<<"nodata value="<<nodata_value<<endl;
    DiskImageView<PixelGray<float> > isis_view( rsrc );
    int width = isis_view.cols();
    int height = isis_view.rows();
    cout<<"width="<<width<<" "<<"height="<<height<<endl;
    camera::IsisCameraModel model(cubFilename);

    InterpolationView<EdgeExtensionView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, ConstantEdgeExtension>, BilinearInterpolation> interpImg
    = interpolate(edge_extend(isis_view,ConstantEdgeExtension()), BilinearInterpolation());
    
    unsigned int index;
    float scaleFactor;
   

    //compute the scale factor between LOLA reflectance and image.
    scaleFactor = ComputeScaleFactor(trackPts);
    cout<<"done computing the scaling factor: "<<scaleFactor<<endl;

    for (index = 0; index < initTransfArray.size(); index++){
      
      errorArray[index] = 0.0;
      numValidPts[index] = 0;
      for (ti = 0; ti < trackPts.size(); ti++){//for each track

	for (si = 0; si < trackPts[ti].size(); si++){//for each shot
          
          //if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0)){
          if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) &&(trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
           
            //float weight = trackPts[ti][si].weightRefl;
        
            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 1/*3*/){//center point of a valid shot
    
                 
		iA = (int) floor(initTransfArray[index][0]*trackPts[ti][si].imgPt[li].x 
                               + initTransfArray[index][1]*trackPts[ti][si].imgPt[li].y 
                               + initTransfArray[index][2]);
		jA = (int) floor(initTransfArray[index][3]*trackPts[ti][si].imgPt[li].x 
                               + initTransfArray[index][4]*trackPts[ti][si].imgPt[li].y 
                               + initTransfArray[index][5]);
                            
                //check (iA,jA) are inside the image!
		if ( ( iA >= 1) && ( iA < height-1) && ( jA >= 1) && ( jA < width-1)){ 
                  //if (interpImg(jA, iA)!=nodata_value){//valid value     
                  if ((float)interpImg(jA, iA)>nodata_value){//valid value  
                    //calculate ii & jj relative to the image center
		    I_e_val = (float)interpImg(jA,iA) - minmax(0) - scaleFactor*trackPts[ti][si].reflectance;
		  }
                  else{
                    I_e_val= 0;
                  }
                  //cout<<"I_e_val="<<I_e_val<<endl;
		  errorArray[index] += abs(I_e_val);
                  numValidPts[index]++;
		}
		
	      }
	    }
	  }
	}
      }
    }

    //determine the best image transform for which the matching error is minimized 
    float minError = errorArray[0];
    int bestIndex = 0;
    finalTransfArray.resize(1);
    matchingErrorArray.resize(1);
    finalTransfArray[0] = initTransfArray[0];
    for (index = 0; index < initTransfArray.size(); index++){
      if (errorArray[index] < minError){
	minError = errorArray[index]/numValidPts[index];
        finalTransfArray[0] = initTransfArray[index];
        bestIndex = index;
      }  
    }

    cout <<"finalTransform: "<<finalTransfArray[0]<<endl;
    cout<<"minError= "<<minError<<endl;
    cout<<"bestIndex="<<bestIndex<<endl; 

    matchingErrorArray[0] = minError;
    
}

//image to lidar coregistration
void UpdateMatchingParamsFromCub(vector<vector<LOLAShot> > &trackPts, string cubFilename,  
			         ModelParams modelParams,  int numMaxIter, 
			         vector<Vector<float, 6> >initTransfArray, vector<float> &initErrorArray, 
                                 vector<Vector<float, 6> >&finalTransfArray, vector<float> &errorArray, Vector2 &centroid)
{

  printf("LIMA matching single processor\n");

  Vector2 minmax = ComputeMinMaxValuesFromCub(cubFilename);

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodata_value = rsrc->nodata_read();
  cout<<"nodata value="<<nodata_value<<endl;
  
  
  DiskImageView<PixelGray<float> > DRG( rsrc );
  int width = DRG.cols();
  int height = DRG.rows();
  cout<<"width="<<width<<" "<<"height="<<height<<endl;

  
  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<float> >   interpDRG = interpolate(edge_extend(DRG,
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing the scale factor...";
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);
  cout<<"scaleFactor="<<scaleFactor<<endl;
  cout<<"done."<<endl;

  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(cubFilename);

  std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_y_deriv.tif";

  cout<<xDerivFilename<<endl;
  cout<<yDerivFilename<<endl;

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  
  float i_C = 0.0;
  float j_C = 0.0;
  int numPts = 0;

  //compute the  i_C and  j_C new way - START
  for (unsigned int ti = 0; ti < trackPts.size(); ti++){
   for (unsigned int si = 0; si < trackPts[ti].size(); si++){
     if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
       for (unsigned int li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot
	 if (trackPts[ti][si].LOLAPt[li].s == 1/*3*/){//center point of a valid shot
            if ((float)interpDRG(trackPts[ti][si].imgPt[li].x, trackPts[ti][si].imgPt[li].y)!=nodata_value){
               i_C = i_C+ trackPts[ti][si].imgPt[li].x;
               j_C = j_C+ trackPts[ti][si].imgPt[li].y;
               numPts++;
	    }
	  }
        }
      }
    }
  }
  i_C = i_C/numPts;
  j_C = j_C/numPts;
  centroid(0) = i_C;
  centroid(1) = j_C;


  int row_max, col_max;
  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
    
  errorArray.resize(initTransfArray.size());
  finalTransfArray.resize(initTransfArray.size());
  vector<int> numValidPts;
  numValidPts.resize(initTransfArray.size());
  
  //copy initTransfArray into finalTransfArray.
  for (unsigned int index = 0; index < initTransfArray.size(); index++){
    for (unsigned int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
    cout <<"InitTransform: "<<initTransfArray[index]<<endl;
  }

  for (unsigned int index = 0; index < initTransfArray.size(); index++){
  
    unsigned int ti,si,li;//indices for tracks, shots and lola points respectively
    int iA, jA;
    int ii, jj;
    int iter;
    float I_x_sqr, I_x_I_y, I_y_sqr; 
    float I_y_val, I_x_val;
    float I_e_val ;
    Matrix<float,6,6> rhs;
    Vector<float,6> lhs;

    iA = 0;
    jA = 0;
    iter = 0;

    while( iter <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs(i_RHS,j_RHS) = 0.0;
        }
        lhs(i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;
      numValidPts[index] = 0;
      
      for (ti = 0; ti < trackPts.size(); ti++){

        for (si = 0; si < trackPts[ti].size(); si++){

          //if ((trackPts[ti][si].valid ==1) && (trackPts[ti][si].reflectance !=/*0*/-1)){
	  if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) && (trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
           
            float weight = trackPts[ti][si].weightRefl;
            //cout<<"weight="<<weight<<endl;        

            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 1/*3*/){//center point of a valid shot
    
		iA = (int) floor(finalTransfArray[index][0]*(trackPts[ti][si].imgPt[li].x-i_C)+ 
                                 finalTransfArray[index][1]*(trackPts[ti][si].imgPt[li].y-j_C)+
                                 finalTransfArray[index][2]+i_C);
		jA = (int) floor(finalTransfArray[index][3]*(trackPts[ti][si].imgPt[li].x-i_C)+ 
                                 finalTransfArray[index][4]*(trackPts[ti][si].imgPt[li].y-j_C)+ 
                                 finalTransfArray[index][5]+j_C);
                
                // check (iA,jA) are inside the image!
		if ( ( iA >= 1) && ( iA < row_max-1) && ( jA >= 1) && ( jA < col_max-1)){ 
	
                  if (((float)interpDRG(jA,iA)!=nodata_value) && (fabs(x_deriv(jA, iA))< 10000) && (fabs(y_deriv(jA, iA))< 10000) ){
 
		    // calculate ii & jj relative to the image center
		    
                    ii = iA - i_C;
		    jj = jA - j_C;
		    
		    I_e_val = (float)interpDRG(jA,iA)-minmax(0) - scaleFactor*trackPts[ti][si].reflectance;
		   
                    /*
                    float robustWeight;
                    float b_sqr = 10;//0.01;
                    if (I_e_val == 0){
		       float tmp = 0.001;
                       robustWeight = sqrt(b_sqr*log(1+(tmp*tmp)/b_sqr))/tmp;
                    }
		    else{
                        robustWeight = sqrt(b_sqr*log(1+(I_e_val*I_e_val)/b_sqr))/fabs(I_e_val);
                    }
                                
                    // We combine the error value with the derivative and
                    // add this to the update equation.
			      	      
		    float weight = robustWeight*trackPts[ti][si].weightRefl;
		    */
                     
                    //cout<<I_e_val<<endl;
                    //start here
                    
		    I_x_val = x_deriv(jA,iA)*weight; 
		    I_y_val = y_deriv(jA,iA)*weight; 
              
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
		  //end here
                    
		  }
                  else{
                    I_e_val= 0;
                  }

		  errorArray[index] += abs(I_e_val);
                  numValidPts[index]++;

		}
	      }           
            }// end of if statement: inside image
          }// end of if statement: valid reflectance  
        }// end of for loop over i
      }//end of loop over k

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

      printRHS(rhs, index, iter);
     
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

      cout<<"numValidPts="<<numValidPts[index]<<endl; 
      errorArray[index] = errorArray[index]/numValidPts[index];
      printLHS_Error(lhs, errorArray[index], index, iter);

      finalTransfArray[index] += lhs;
  
      iter ++;
 
    }
   
  }//index loop ends here


  //TODO: check to make sure the error is smaller than after the initialization stage.
} 


































//determines the best finalTransfArray from all initTransfArrays
//it assumes that each image is transformed by one affine transform
//in the future we can investigate the use of one affine transform per track. 
void InitMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
			ModelParams modelParams, CoregistrationParams coregistrationParams,  
			vector<Vector<float, 6> >initTransfArray, Vector<float, 6> &finalTransf, 
			float *matchingError)
{
    unsigned int ti, si, li;
    int jA, iA;
    float I_e_val;
  
    vector<float>errorArray;
    errorArray.resize(initTransfArray.size());
    DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
    GeoReference DRGGeo;
    read_georeference(DRGGeo, DRGFilename);

    cout <<"Interpolating the image ...";
    ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
									 ConstantEdgeExtension()),
							     BilinearInterpolation());
    cout<<"done."<<endl;

    unsigned int index;
    float scaleFactor;
    int row_max = DRG.rows();
    int col_max = DRG.cols();
  
    //compute the scale factor between LOLA reflectance and image.
    scaleFactor = ComputeScaleFactor(trackPts);
    cout<<"done computing the scaling factor"<<endl;

    for (index = 0; index < initTransfArray.size(); index++){
      
      errorArray[index] = 0.0;
      
      for (ti = 0; ti < trackPts.size(); ti++){//for each track

	for (si = 0; si < trackPts[ti].size(); si++){//for each shot
          
          //cout<<ti<<" "<<si<<endl;

          if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != /*0*/-1)){
           
            // float weight = trackPts[ti][si].weightRefl;
        
            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 3){//center point of a valid shot
    
                 
		iA = (int) floor(initTransfArray[index][0]*trackPts[ti][si].imgPt[li].x 
                               + initTransfArray[index][1]*trackPts[ti][si].imgPt[li].y 
                               + initTransfArray[index][2]);
		jA = (int) floor(initTransfArray[index][3]*trackPts[ti][si].imgPt[li].x 
                               + initTransfArray[index][4]*trackPts[ti][si].imgPt[li].y 
                               + initTransfArray[index][5]);
                
                //check (iA,jA) are inside the image!
		if ( ( iA >= 0) && ( iA < row_max) && ( jA >= 0) && ( jA < col_max)){ 
		  if (is_valid(interpDRG(jA, iA)) ){
                    //calculate ii & jj relative to the image center
		    I_e_val = interpDRG(jA,iA) - scaleFactor*trackPts[ti][si].reflectance;
		  }
                  else{
                    I_e_val= 0;
                  }

		  errorArray[index] += abs(I_e_val);
		}
		
	      }
	    }
	  }
	}
      }
    }

    //determine the best image transform for which the matching error is minimized 
    float minError = errorArray[0];
    finalTransf = initTransfArray[0];
    for (index = 0; index < initTransfArray.size(); index++){
      if (errorArray[index] < minError){
	minError = errorArray[index];
        finalTransf = initTransfArray[index];
      }  
    } 
    *matchingError = minError;
}



//image to lidar coregistration
void UpdateMatchingParams(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
			  ModelParams modelParams,  int numMaxIter, 
			  vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
			  vector<float> &errorArray )
{

  printf("LIMA matching single processor\n");

  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing the scale factor...";
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);
  cout<<"done."<<endl;

  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DRGFilename);

  std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_y_deriv.tif";

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  

  //float xx,yy;
  int row_max, col_max;

  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
   
  cout << "image size: rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C = row_max/2;
  int j_C = col_max/2;

  //working with the entire CVS file significant errors where seen at index 126 - set index to this by hand to debugg this test case
 /*
  1. we see huge problems from the get go - massive lhs numbers, where does these come from in comparison
  2. need to explore & instrument the code further
  3. we have some stable/unstable pairs to explore (106 stable, 109 not) & ( 112 stable, 115 not)
 
  A couple of thoughts - we rezero lhs & rhs every iteration this should not be the problem.  Could the problem simple be a scaling issue where the size of the error terms on the right & left hand side - caused by the great number of points - causes numerical stability issues for solving the matrix?
*/

  //copy initTransfArray into finalTransfArray.
  for (unsigned int index = 0; index < initTransfArray.size(); index++){
    for (unsigned int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
  }

  for (unsigned int index = 0; index < initTransfArray.size(); index++){

    cout << "index = "<< index << endl;
  
  
    unsigned int ti,si,li;//indices for tracks, shots and lola points respectively
    int iA, jA;
    int ii, jj;
    int iter;
    float I_x_sqr, I_x_I_y, I_y_sqr; 
    float I_y_val, I_x_val;
    float I_e_val ;
    Matrix<float,6,6> rhs;
    Vector<float,6> lhs;
 
    iA = 0;
    jA = 0;
    iter = 0;

    while( iter <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs(i_RHS,j_RHS) = 0.0;
        }
        lhs(i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;
      
      for (ti = 0; ti < trackPts.size(); ti++){

        for (si = 0; si < trackPts[ti].size(); si++){

          //if ((trackPts[ti][si].valid ==1) && (trackPts[ti][si].reflectance !=0)){
          if ((trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0) &&(trackPts[ti][si].reflectance != -1)){//valid track and non-zero reflectance
           
            float weight = trackPts[ti][si].weightRefl;
        
            for (li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot

              if (trackPts[ti][si].LOLAPt[li].s == 3){//center point of a valid shot
    
		iA = (int) floor(finalTransfArray[index][0]*trackPts[ti][si].imgPt[li].x + 
                                 finalTransfArray[index][1]*trackPts[ti][si].imgPt[li].y + 
                                 finalTransfArray[index][2]);
		jA = (int) floor(finalTransfArray[index][3]*trackPts[ti][si].imgPt[li].x + 
                                 finalTransfArray[index][4]*trackPts[ti][si].imgPt[li].y + 
                                 finalTransfArray[index][5]);
                
                // check (iA,jA) are inside the image!
		if ( ( iA >= 0) && ( iA < row_max) && ( jA >= 0) && ( jA < col_max)){ 
		  if (is_valid(interpDRG(jA, iA)) ){
                  // calculate ii & jj relative to the image center
		    ii = iA - i_C;
		    jj = jA - j_C;
		    I_e_val = interpDRG(jA,iA) - scaleFactor*trackPts[ti][si].reflectance;
		  }
                  else{
                    I_e_val= 0;
                  }

		  errorArray[index] += abs(I_e_val);
		  
                
		  /*
                  float robustWeight;
                  float b_sqr = 0.0001;
                  if (I_e_val == 0){
		     float tmp = 0.001;
                     robustWeight = sqrt(b_sqr*log(1+(tmp*tmp)/b_sqr))/tmp;
                  }
		  else{
                     robustWeight = sqrt(b_sqr*log(1+(I_e_val*I_e_val)/b_sqr))/fabs(I_e_val);
                  }
                                
                  // We combine the error value with the derivative and
                  // add this to the update equation.
			      	      
		  float weight = spatialWeights(ii+kern_half_width, jj+kern_half_height)*robustWeight;
            
		  */


		  //calculate numerical dirivatives (ii,jj).
                
		  I_x_val = x_deriv(jA,iA)*weight; 
		  I_y_val = y_deriv(jA,iA)*weight; 

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

		}
	      }           
            }// end of if statement: inside image
          }// end of if statement: valid reflectance  
        }// end of for loop over i
      }//end of loop over k

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

   
      printRHS(rhs, index, iter);
     
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


      printLHS_Error(lhs, errorArray[index], index, iter);

      finalTransfArray[index] += lhs;
  
      iter ++;
 
    }
   
  }//index loop ends here

} 

#if 0
//this is the MP-multi processor version of the above function
void UpdateMatchingParamsLIMA_MP(vector<vector<LOLAShot> > &trackPts, string DRGFilename,  
				 ModelParams modelParams, CoregistrationParams coregistrationParams,  
				 vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
				 vector<float> &errorArray )
{

  int numMaxIter = coregistrationParams.maxNumIter;

  DiskImageView<PixelGray<uint8> >   DRG(DRGFilename);
  GeoReference DRGGeo;
  read_georeference(DRGGeo, DRGFilename);

  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<uint8> >   interpDRG = interpolate(edge_extend(DRG.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing the scale factor...";
  float scaleFactor;
  scaleFactor = ComputeScaleFactor(trackPts);
  cout<<"done."<<endl;

  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DRGFilename);

  std::string xDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_less3_from_filename(temp) + "_y_deriv.tif";

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DRG),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  

  int row_max, col_max;

  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
   
  cout << "image size: rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C = row_max/2;
  int j_C = col_max/2;


  //copy initTransfArray into finalTransfArray.
  for (int index = 0; index < initTransfArray.size(); index++){
    for (int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
  }

   vector<int> ti;
   vector<int> si;
   vector<int> li;
   vector<int> iter;
   vector<int> ii;
   vector<int> jj;
   vector<int> iA;
   vector<int> jA;
   vector<float> I_e_val;
   vector<float> I_x_val;
   vector<float> I_y_val;
   vector<float> I_x_sqr;
   vector<float> I_x_I_y;
   vector<float> I_y_sqr;
   vector<Matrix<float, 6, 6> > rhs;
   vector<Vector<float,6> > lhs; 

   ti.resize(initTransfArray.size());
   si.resize(initTransfArray.size());
   li.resize(initTransfArray.size());
   iter.resize(initTransfArray.size());
   ii.resize(initTransfArray.size());
   jj.resize(initTransfArray.size());
   iA.resize(initTransfArray.size());
   jA.resize(initTransfArray.size());
   I_e_val.resize(initTransfArray.size());
   I_x_val.resize(initTransfArray.size());
   I_y_val.resize(initTransfArray.size());
   I_x_sqr.resize(initTransfArray.size());
   I_x_I_y.resize(initTransfArray.size());
   I_y_sqr.resize(initTransfArray.size());
   rhs.resize(initTransfArray.size());
   lhs.resize(initTransfArray.size());


   int nthreads, tid, chunk; 
   chunk = CHUNKSIZE;
   for (int index = 0; index< initTransfArray.size(); index++){
        iA[index] = 0;
        jA[index] = 0;
        iter[index] = 0;
   }

#pragma omp parallel shared(ti,si, li, iter, ii, jj, iA, jA, I_e_val, I_x_val, I_y_val, I_x_sqr, I_x_I_y, I_y_sqr, rhs, lhs, nthreads,chunk) private(index,tid)
{
  
  tid = omp_get_thread_num();
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
  }
  printf("Thread %d starting...\n",tid);


#pragma omp for schedule(dynamic,chunk)

  for (int index = 0; index < initTransfArray.size(); index++){

    cout << "index = "<< index << endl;

    while( iter[index] <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs[index](i_RHS,j_RHS) = 0.0;
        }
        lhs[index](i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;

      for (ti[index] = 0; ti[index] < trackPts.size(); ti[index]++){

        for (si[index] = 0; si[index] < trackPts[ti[index]].size(); si[index]++){

          if ((trackPts[ti[index]][si[index]].valid ==1) && (trackPts[ti[index]][si[index]].reflectance !=0)){
           
            float weight = trackPts[ti[index]][si[index]].weightRefl;
        
            for (li[index] = 0; li[index] < trackPts[ti[index]][si[index]].LOLAPt.size(); li[index]++){//for each point of a LOLA shot

              if (trackPts[ti[index]][si[index]].LOLAPt[li[index]].s == 3){//center point of a valid shot
    
		iA[index] = (int) floor(finalTransfArray[index][0]*trackPts[ti[index]][si[index]].imgPt[li[index]].x 
                          + finalTransfArray[index][1]*trackPts[ti[index]][si[index]].imgPt[li[index]].y + finalTransfArray[index][2]);
		jA[index] = (int) floor(finalTransfArray[index][3]*trackPts[ti[index]][si[index]].imgPt[li[index]].x 
                          + finalTransfArray[index][4]*trackPts[ti[index]][si[index]].imgPt[li[index]].y + finalTransfArray[index][5]);
                

                // check (iA,jA) are inside the derivative image!
		if ( ( iA[index] >= 0) && ( iA[index] < row_max) && ( jA[index] >= 0) && ( jA[index] < col_max)){

                    	  
                  // calculate ii & jj relative to the image center
		  ii[index] = iA[index] - i_C;
		  jj[index] = jA[index] - j_C;

		  I_e_val[index] = interpDRG(jA[index],iA[index]) - scaleFactor*trackPts[ti[index]][si[index]].reflectance;
		  errorArray[index] += abs(I_e_val[index]);

		  //calculate numerical derivatives (ii,jj).
                  //printf("weight = %f\n", weight);

		  I_x_val[index] = x_deriv(jA[index],iA[index])*weight; 
		  I_y_val[index] = y_deriv(jA[index],iA[index])*weight; 

		  I_x_I_y[index] = I_x_val[index]*I_y_val[index];        
		  I_x_sqr[index] = I_x_val[index]*I_x_val[index];
		  I_y_sqr[index] = I_y_val[index]*I_y_val[index];
                  
                  //printf("I_x_val=%f, I_y_val=%f, I_e_val=%f\n", I_x_val[index],  I_y_val[index], I_e_val[index]);     
		  // Left hand side
		  lhs[index](0) += ii[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](1) += jj[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](2) +=                              I_x_val[index] * I_e_val[index];
		  lhs[index](3) += ii[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](4) += jj[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](5) +=                              I_y_val[index] * I_e_val[index];

		  // Right Hand Side UL
		  rhs[index](0,0) += ii[index]*ii[index] * I_x_sqr[index];
		  rhs[index](0,1) += ii[index]*jj[index] * I_x_sqr[index];
		  rhs[index](0,2) += ii[index]           * I_x_sqr[index];
		  rhs[index](1,1) += jj[index]*jj[index] * I_x_sqr[index];
		  rhs[index](1,2) += jj[index]           * I_x_sqr[index];
		  rhs[index](2,2) +=                       I_x_sqr[index];

		  // Right Hand Side UR
		  rhs[index](0,3) += ii[index]*ii[index] * I_x_I_y[index];
		  rhs[index](0,4) += ii[index]*jj[index] * I_x_I_y[index];
		  rhs[index](0,5) += ii[index]           * I_x_I_y[index];
		  rhs[index](1,4) += jj[index]*jj[index] * I_x_I_y[index];
		  rhs[index](1,5) += jj[index]           * I_x_I_y[index];
		  rhs[index](2,5) +=                       I_x_I_y[index];

		  // Right Hand Side LR
		  rhs[index](3,3) += ii[index]*ii[index] * I_y_sqr[index];
		  rhs[index](3,4) += ii[index]*jj[index] * I_y_sqr[index];
		  rhs[index](3,5) += ii[index]           * I_y_sqr[index];
		  rhs[index](4,4) += jj[index]*jj[index] * I_y_sqr[index];
		  rhs[index](4,5) += jj[index]           * I_y_sqr[index];
		  rhs[index](5,5) +=                       I_y_sqr[index];

		}
	      }           
            }// end of if statement: inside image
          }// end of if statement: valid reflectance  
        }// end of for loop over i
      }//end of loop over k

      // Fill in symmetric entries
      rhs[index](1,0) = rhs[index](0,1);
      rhs[index](2,0) = rhs[index](0,2);
      rhs[index](2,1) = rhs[index](1,2);
      rhs[index](1,3) = rhs[index](0,4);
      rhs[index](2,3) = rhs[index](0,5);
      rhs[index](2,4) = rhs[index](1,5);
      rhs[index](3,0) = rhs[index](0,3);
      rhs[index](3,1) = rhs[index](1,3);
      rhs[index](3,2) = rhs[index](2,3);
      rhs[index](4,0) = rhs[index](0,4);
      rhs[index](4,1) = rhs[index](1,4);
      rhs[index](4,2) = rhs[index](2,4);
      rhs[index](4,3) = rhs[index](3,4);
      rhs[index](5,0) = rhs[index](0,5);
      rhs[index](5,1) = rhs[index](1,5);
      rhs[index](5,2) = rhs[index](2,5);
      rhs[index](5,3) = rhs[index](3,5);
      rhs[index](5,4) = rhs[index](4,5);

      //printf ("before\n");
      //printRHS(rhs[index], index, iter[index]);
      //printLHS_Error(lhs[index], errorArray[index], index, iter[index]);

      try {
        solve_symmetric_nocopy(rhs[index],lhs[index]);
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


      printLHS_Error(lhs[index], errorArray[index], index, iter[index]);

      finalTransfArray[index] += lhs[index];
  
      iter[index] ++;
 
    }
   
  }//index loop ends here

 }//openMP


} 
#endif


//this is the MP-multi processor version of the above function
void UpdateMatchingParamsLIDEM_MP(vector<vector<LOLAShot> > &trackPts, string DEMFilename,  
				  ModelParams modelParams,  CoregistrationParams coregistrationParams, 
				  vector<Vector<float, 6> >initTransfArray, vector<Vector<float, 6> >&finalTransfArray, 
				  vector<float> &errorArray )
{
#if 0
  int numMaxIter = coregistrationParams.maxNumIter;

  DiskImageView<PixelGray<float> >   DEM(DEMFilename);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, DEMFilename);

  cout <<"Interpolating the image ...";
  ImageViewRef<PixelGray<float> >   interpDEM = interpolate(edge_extend(DEM.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());
  cout<<"done."<<endl;


  cout << "Computing/Reading the derivative image..."; 
  //to be used in the RELEASE!- DO NOT remove it!
  //DiskCacheImageView<float> x_deriv = derivative_filter(DRG,1,0);
  //DiskCacheImageView<float> y_deriv = derivative_filter(DRG,0,1);

  std::string temp = sufix_from_filename(DEMFilename);

  std::string xDerivFilename = "../results" + prefix_from_filename(temp) + "_x_deriv.tif";
  std::string yDerivFilename = "../results" + prefix_from_filename(temp) + "_y_deriv.tif";

  if ( !boost::filesystem::exists( xDerivFilename ) ) {
    cout << "Computing the x_derivative ..." << endl;
    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DEM),1,0);
    DiskImageResourceGDAL rsrc( xDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the x_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
    cin.get();
    cout<<"done."<<endl;
  }
  cout << "Reading in the x_deriv from "<<xDerivFilename<<" ..." << endl;

  DiskImageView<float> x_deriv( xDerivFilename );

  if ( !boost::filesystem::exists( yDerivFilename ) ) {
     cout << "Computing the y_derivative ..." << endl;

    ImageViewRef<float> temp = derivative_filter(pixel_cast<float>(DEM),0,1);
    DiskImageResourceGDAL rsrc( yDerivFilename,
        temp.format(), Vector2i(512,512) );
    cout << "Writing the y_derivative" << endl;
    block_write_image(rsrc, temp,
        TerminalProgressCallback("asp", "Derivative:") );
  }
 
  cout << "Reading in the y_deriv from "<<yDerivFilename<<" ... "<< endl;
  DiskImageView<float> y_deriv( yDerivFilename );
  cout<<"done."<<endl;
  

  int row_max, col_max;

  row_max = x_deriv.rows();
  col_max = x_deriv.cols();
   
  cout << "image size: rows: " << row_max << ", cols:"<< col_max << endl;
  int i_C = row_max/2;
  int j_C = col_max/2;


  //copy initTransfArray into finalTransfArray.
  for (int index = 0; index < initTransfArray.size(); index++){
    for (int i = 0; i < 6; i++){
      finalTransfArray[index][i]=initTransfArray[index][i];
    }
  }

  
   vector<int> ti;
   vector<int> si;
   vector<int> li;
   vector<int> iter;
   vector<int> ii;
   vector<int> jj;
   vector<int> iA;
   vector<int> jA;
   vector<float> I_e_val;
   vector<float> I_x_val;
   vector<float> I_y_val;
   vector<float> I_x_sqr;
   vector<float> I_x_I_y;
   vector<float> I_y_sqr;
   vector<Matrix<float, 6, 6> > rhs;
   vector<Vector<float,6> > lhs; 

   ti.resize(initTransfArray.size());
   si.resize(initTransfArray.size());
   li.resize(initTransfArray.size());
   iter.resize(initTransfArray.size());
   ii.resize(initTransfArray.size());
   jj.resize(initTransfArray.size());
   iA.resize(initTransfArray.size());
   jA.resize(initTransfArray.size());
   I_e_val.resize(initTransfArray.size());
   I_x_val.resize(initTransfArray.size());
   I_y_val.resize(initTransfArray.size());
   I_x_sqr.resize(initTransfArray.size());
   I_x_I_y.resize(initTransfArray.size());
   I_y_sqr.resize(initTransfArray.size());
   rhs.resize(initTransfArray.size());
   lhs.resize(initTransfArray.size());

   int index;
   int nthreads, tid, chunk;
     
   chunk = CHUNKSIZE;
   for (int index = 0; index< initTransfArray.size(); index++){
        iA[index] = 0;
        jA[index] = 0;
        iter[index] = 0;
   }

#pragma omp parallel shared(ti, si, li, iter, ii, jj, iA, jA, I_e_val, I_x_val, I_y_val, I_x_sqr, I_x_I_y, I_y_sqr, rhs, lhs, nthreads,chunk) private(index,tid)
{
  
  tid = omp_get_thread_num();
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
  }
  printf("Thread %d starting...\n",tid);


#pragma omp for schedule(dynamic,chunk)

  for (int index = 0; index < initTransfArray.size(); index++){
    //for (index = 124; index < 128; index++){
 
    while( iter[index] <= numMaxIter){ //gradient descent => optimal transform

      //reset rhs & lhs
      for(int i_RHS = 0; i_RHS < 6; i_RHS++){
        for(int j_RHS = 0; j_RHS < 6; j_RHS++){
          rhs[index](i_RHS,j_RHS) = 0.0;
        }
        lhs[index](i_RHS) = 0.0;
      }
 
      //reset the error;
      errorArray[index] = 0.0;

      for (ti[index] = 0; ti[index] < trackPts.size(); ti[index]++){

        for (si[index] = 0; si[index] < trackPts[ti[index]].size(); si[index]++){

          if (trackPts[ti[index]][si[index]].valid ==1){
           
            //weight = trackPts[k][i].weight_prd;
        
            for (li[index] = 0; li[index] < trackPts[ti[index]][si[index]].LOLAPt.size(); li[index]++){//for each point of a LOLA shot

     
		iA[index] = (int) floor(finalTransfArray[index][0]*trackPts[ti[index]][si[index]].DEMPt[li[index]].x 
                          + finalTransfArray[index][1]*trackPts[ti[index]][si[index]].DEMPt[li[index]].y + finalTransfArray[index][2]);
		jA[index] = (int) floor(finalTransfArray[index][3]*trackPts[ti[index]][si[index]].DEMPt[li[index]].x 
                          + finalTransfArray[index][4]*trackPts[ti[index]][si[index]].DEMPt[li[index]].y + finalTransfArray[index][5]);
                
                // check (iA,jA) are inside the image!
		if ( ( iA[index] >= 0) && ( iA[index] < row_max) && ( jA[index] >= 0) && ( jA[index] < col_max)){
        	  
                  // calculate ii & jj relative to the image center
		  ii[index] = iA[index] - i_C;
		  jj[index] = jA[index] - j_C;
                 
		  I_e_val[index] = interpDEM(jA[index],iA[index]) - trackPts[ti[index]][si[index]].LOLAPt[li[index]].coords[2];
               
		  errorArray[index] += abs(I_e_val[index]);

		  //calculate numerical dirivatives (ii,jj).
		  I_x_val[index] = x_deriv(jA[index],iA[index]); 
		  I_y_val[index] = y_deriv(jA[index],iA[index]); 

		  I_x_I_y[index] = I_x_val[index]*I_y_val[index];        
		  I_x_sqr[index] = I_x_val[index]*I_x_val[index];
		  I_y_sqr[index] = I_y_val[index]*I_y_val[index];

		  // Left hand side
		  lhs[index](0) += ii[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](1) += jj[index] * I_x_val[index] * I_e_val[index];
		  lhs[index](2) +=                              I_x_val[index] * I_e_val[index];
		  lhs[index](3) += ii[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](4) += jj[index] * I_y_val[index] * I_e_val[index];
		  lhs[index](5) +=                              I_y_val[index] * I_e_val[index];

		  // Right Hand Side UL
		  rhs[index](0,0) += ii[index]*ii[index] * I_x_sqr[index];
		  rhs[index](0,1) += ii[index]*jj[index] * I_x_sqr[index];
		  rhs[index](0,2) += ii[index]           * I_x_sqr[index];
		  rhs[index](1,1) += jj[index]*jj[index] * I_x_sqr[index];
		  rhs[index](1,2) += jj[index]           * I_x_sqr[index];
		  rhs[index](2,2) +=                       I_x_sqr[index];

		  // Right Hand Side UR
		  rhs[index](0,3) += ii[index]*ii[index] * I_x_I_y[index];
		  rhs[index](0,4) += ii[index]*jj[index] * I_x_I_y[index];
		  rhs[index](0,5) += ii[index]           * I_x_I_y[index];
		  rhs[index](1,4) += jj[index]*jj[index] * I_x_I_y[index];
		  rhs[index](1,5) += jj[index]           * I_x_I_y[index];
		  rhs[index](2,5) +=                       I_x_I_y[index];

		  // Right Hand Side LR
		  rhs[index](3,3) += ii[index]*ii[index] * I_y_sqr[index];
		  rhs[index](3,4) += ii[index]*jj[index] * I_y_sqr[index];
		  rhs[index](3,5) += ii[index]           * I_y_sqr[index];
		  rhs[index](4,4) += jj[index]*jj[index] * I_y_sqr[index];
		  rhs[index](4,5) += jj[index]           * I_y_sqr[index];
		  rhs[index](5,5) +=                       I_y_sqr[index];

		}//end of if statemenet: check (iA,jA) are inside the image! 
            }// end of for loop over li
          }// end of if statement: valid track and reflectance  
        }// end of for loop over si
      }//end of loop over ti

      // Fill in symmetric entries
      rhs[index](1,0) = rhs[index](0,1);
      rhs[index](2,0) = rhs[index](0,2);
      rhs[index](2,1) = rhs[index](1,2);
      rhs[index](1,3) = rhs[index](0,4);
      rhs[index](2,3) = rhs[index](0,5);
      rhs[index](2,4) = rhs[index](1,5);
      rhs[index](3,0) = rhs[index](0,3);
      rhs[index](3,1) = rhs[index](1,3);
      rhs[index](3,2) = rhs[index](2,3);
      rhs[index](4,0) = rhs[index](0,4);
      rhs[index](4,1) = rhs[index](1,4);
      rhs[index](4,2) = rhs[index](2,4);
      rhs[index](4,3) = rhs[index](3,4);
      rhs[index](5,0) = rhs[index](0,5);
      rhs[index](5,1) = rhs[index](1,5);
      rhs[index](5,2) = rhs[index](2,5);
      rhs[index](5,3) = rhs[index](3,5);
      rhs[index](5,4) = rhs[index](4,5);

   
      //printRHS(rhs[index], index, iter[index]);
     
      try {
        solve_symmetric_nocopy(rhs[index],lhs[index]);
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


      //printLHS_Error(lhs[index], errorArray[index], index, iter[index]);

      finalTransfArray[index] += lhs[index];
  
      iter[index] ++;
 
    }
   
  }//index loop ends here

 }//openMP
#endif

} 


/*
//saves the 3D and image points in a local control network
//will be renamed to save control network
void SaveImagePts(vector<vector<LOLAShot> > &trackPts, Vector<float, 6> finalTransfArray, float error,  string matchResultsFilename)
{

   printf("saving the control network ...\n");
   FILE * fp = fopen(matchResultsFilename.c_str(),"w");

   for (int ti = 0; ti < trackPts.size(); ti++){

      for (int si = 0; si < trackPts[ti].size(); si++){
	if ((trackPts[ti][si].featurePtRefl == 1.0) && (trackPts[ti][si].valid == 1) && (trackPts[ti][si].reflectance != 0)) { 
           for (int li = 0; li < trackPts[ti][si].LOLAPt.size(); li++){//for each point of a LOLA shot
              if (trackPts[ti][si].LOLAPt[li].s == 3){//center point of a valid shot
		int i = (int) floor(finalTransfArray[0]*trackPts[ti][si].imgPt[li].x + finalTransfArray[1]*trackPts[ti][si].imgPt[li].y + finalTransfArray[2]);
		int j = (int) floor(finalTransfArray[3]*trackPts[ti][si].imgPt[li].x + finalTransfArray[4]*trackPts[ti][si].imgPt[li].y + finalTransfArray[5]);
		fprintf(fp,"x: %f, y: %f, z: %f, i: %d,  j: %d\n", 
                        trackPts[ti][si].LOLAPt[li].coords[0], 
			trackPts[ti][si].LOLAPt[li].coords[1], 
                        trackPts[ti][si].LOLAPt[li].coords[2], i, j);
	      } 
	   }
        } 
      }
   }

   fclose(fp);
   printf("done.\n");
}
*/
