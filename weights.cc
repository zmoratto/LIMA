#include "weights.h"

//determine the weights of a LOLA track
//it should be moved in weights.cc
void ComputeWeights( vector< vector<LOLAShot> > & trackPts, int &halfWindow, float & take_p, int&  numValid, float & take_thresh,string filename ){
  InterpolateInvalidPoint( trackPts);
  FindValidPoints( trackPts, halfWindow,  numValid);
  ComputeGradient( trackPts, halfWindow);
  ComputeSalientFeatures( trackPts,  take_p, numValid, take_thresh);
  SaveWeights( trackPts, filename);
}


//this function loops through each set of track points separatly, 
//if we see the following pattern of points (valid, in-valid, valid)
//the non-valid point will have it's reflectance set to the average of the two valid points
//
//above non-valid was defined as reflectance = -1 
int InterpolateInvalidPoint(vector<vector< LOLAShot> > & trackPts){

  for(int ti = 0; ti < trackPts.size() ; ti++ ){ //index track

    //debugg inner loop - correct invalid points
    for(int si = 1; si < trackPts[ti].size()-1; si++ ){
    
      //printf("reflectance at %d = %f\n",iPt,track_pts[iTrk][iPt].reflectance);
      if(trackPts[ti][si].reflectance < 0){
        if ((trackPts[ti][si-1].reflectance > 0 ) && (trackPts[ti][si+1].reflectance > 0)){
           trackPts[ti][si].reflectance = (trackPts[ti][si-1].reflectance + trackPts[ti][si+1].reflectance )/2;  
        }
      }
    }
  }

  //printf("\n\ncorrect_invalid_linear = done\n\n");
  return 0;
}


//function loops through points setting the LOLAShot integer field .calc_acp to "1" 
//if the point is valid & there are sufficent valid points on either side of that point compute a gradient.  
//otherwise .calc_acp is set to 0 if either the point is invalid or there are insufficient 
//valid points bordering that point to compute a gradient
//halfWindow = (length of the derivative filter-1)/2
//numValid = the number of points for which the derivative filter is computed
 
int FindValidPoints(vector< vector<LOLAShot > > & trackPts,int halfWindow, int &numValid){


  int t_inv = -1; //trailing invalid
  int l_inv = -1; //leading invalid
  numValid = 0;
  int ti; //track index
  int si; //shot  index
  for (ti = 0; ti < trackPts.size(); ti ++){
  
    int numShots = trackPts[ti].size();
    t_inv = -1;
    l_inv = -1;
  
    for (si = 0; si < numShots; si++){

      //don't compute the derivatives for points close to the track boundary
      if (si >= numShots-halfWindow){
         trackPts[ti][si].calc_acp = false;
         continue;
      }

      if (trackPts[ti][si].reflectance < 0){  
          t_inv = si;
      }
      if (trackPts[ti][si+halfWindow].reflectance < 0){
          l_inv = si + halfWindow;
      }

      if ((si > (halfWindow-1)) && (si < (numShots-halfWindow+1)) ){
        
        if ( ((si-t_inv) > halfWindow) && (l_inv < si)){
           trackPts[ti][si].calc_acp = 1;
           numValid ++;
        }
        else{
           trackPts[ti][si].calc_acp = 0;
        }

      }

      else{
          trackPts[ti][si].calc_acp = 0;
      } 
    
    }
  }
  printf("\n\nEnd FindValidPts: calculations of acceptable points \n\n"); 
  return 0;
}

//calculates the gradient of reflectance using a filter of shape: [ -1x--edge-- 0  +1x--edge--]
int ComputeGradient(vector< vector<LOLAShot > > & trackPts, int halfWindow){
  float sum = 0; // filter response


  cout << "Starting grad filter application\n" << endl;
  

  //build filter -START
  vector<float> f;  
  int windowSize = 2*halfWindow + 1;
  f.resize(windowSize);

  for (int i = 0; i < windowSize ;i++ ){
    if (i < halfWindow){
      f[i] = -1.0; 
    }
    else if(i==halfWindow){
       f[i] = 0.0;
    }
    else{
      f[i] = 1.0;
    }
  }
  //build filter -END

  //apply filter

  //cout << "Get to start of filter application?" << endl;;
  for(int ti = 0; ti < trackPts.size(); ti++){

    int numShots = trackPts[ti].size();

    printf("track_pts[%d].size() = %d\n",ti,numShots);

    for(int si = halfWindow; si < (numShots - halfWindow); si++){

      if( trackPts[ti][si].calc_acp){
        sum = 0;
        for(int place = 0; place < windowSize; place ++){
          sum += f[place]*trackPts[ti][si-halfWindow+place].reflectance;
        }
        trackPts[ti][si].filter_response = sum;
         //printf("track_pts[%d][%d] = %0.1f\n", iTrk , ip, track_pts[iTrk][ip].filter_response);
      }
    }
  }

  //cout << "End of filter application" << endl;
  return 0;
}


// we'll start with [0,1] based on top 5% 
// then allow room for growth
//3.compute the salient features of the LOLA tracks


int ComputeSalientFeatures( vector<vector< LOLAShot > > & trackPts, float topPercent, int& numValid, float& take_thresh){
  // builds a list of the abs(gradient response), 
  // finds the response value that corresponds with the top "take_p" %
  // assignes LOLAShot.weight_lsq = 1 if the response of that point is in the top take_p of all points
  // 
  // 
  //make weights - start with top 'take_p' accross all tracks
  //initial implementation will be slow

  // 0. make a list
  vector<float> list_all_response;
  int numRead = 0;
  
  list_all_response.resize(numValid);

  for(int ti = 0; ti < trackPts.size(); ti++ ){
    //printf("list_all_response.size() = %d\n", list_all_response.size() );

    for( int si = 0; si < trackPts[ti].size(); si++ ){
 
      if(trackPts[ti][si].calc_acp){
        list_all_response[numRead] = abs(trackPts[ti][si].filter_response);
        numRead ++;
      }
    }
  }
 
  // 1. sort list
  sort(list_all_response.begin(),list_all_response.end() );
 
 
  // 2. define the threshold
  int take_point = 0;
  take_point = (int)ceil( (1-topPercent) * numValid);
  take_thresh = 0.0;
  take_thresh = list_all_response[take_point];
  printf("list_all_response.size() = %d\n",(int)list_all_response.size());
  printf("make_weights report: take_point = %d, take_thresh = %f\n",take_point, take_thresh);
 

  // 3. assign values - 
  int count_greater = 0;
  int others = 0;
  for(int ti = 0; ti < trackPts.size(); ti++ ){
    for(int si = 0; si < trackPts[ti].size(); si ++ ){
      //cout << "itr_all = " << itr_all << endl;
      if( abs(trackPts[ti][si].filter_response) >= take_thresh ){
        trackPts[ti][si].weight_lsq = 1.0;  
        count_greater ++;
      }
      else{
        trackPts[ti][si].weight_lsq = 0.05;
        others ++;
      }
    }
  }
  printf("count_greater = %d, others = %d", count_greater, others);
  printf("Make weights function should be complete - take_point = %d, take_thresh = %f\n", take_point, take_thresh);
  return 0;
}

//assign linear weights (shapped like a pyramid) to all points around one of the interest point selected by make_weights
int MakeLinearWeights( vector< vector< LOLAShot > > & trackPts, const int &halfWindow, const int & windowSize){

  float windowSize_f = windowSize;
  float halfWindow_f = halfWindow;

  float to_add = 0.0;
  printf("windowSize_f = %f, halfWindow_f = %f\n", windowSize_f, halfWindow_f);

  for(int ti = 0; ti < trackPts.size(); ti++ ){
    //zero the track before begining 
    for(int si = 0; si < trackPts[ti].size(); si++){
      trackPts[ti][si].weight_prd = 0; // just to begin
    }

    for( int si = 0; si < trackPts[ti].size() ; si++){
      if( trackPts[ti][si].weight_lsq == 1 ){

        //write out using edge - calculate normalization latter
        to_add = 0.0;
        for( int m = 0; m < windowSize; m++){
          to_add = (windowSize - abs(m-halfWindow_f)) / windowSize_f;
          trackPts[ti][si-halfWindow+m].weight_prd += to_add;
        }
      }
    }
  }

  printf("Done computing the linear weights\n");
  return 0;
}

//write out all weight specific variables to a txt file specified by "f_name".  
//writen out is: reflectance, calc_acp, filter_response, weight_lsq, weight_prd
int SaveWeights( vector< vector< LOLAShot> > & trackPts, string & filename){
  
 cout << "Begining the end - write out commensing " << endl;
 FILE* sFile;
 sFile = fopen(filename.c_str(), "w");
 fprintf(sFile, "Write weights:\n");
 for (int ti = 0; ti < trackPts.size(); ti++ ){
   for (int si = 0; si < trackPts[ti].size(); si++ ){
     fprintf( sFile, "ti=  %d si= %d reflectance= %f calc_acp= %d filter_response= %f weight= %f weight_prd= %f \n", 
                     ti, si, trackPts[ti][si].reflectance, trackPts[ti][si].calc_acp, 
                     trackPts[ti][si].filter_response, trackPts[ti][si].weight_lsq, trackPts[ti][si].weight_prd );
   }
 }

 fclose(sFile);
 cout << "Write out completed"<< endl;
}

