#include "weights.h"

//computes and saves the weights of a LOLA track
void ComputeWeights( vector< vector<LOLAShot> >& trackPts, int halfWindow, float topPercent, string filename ){
 
  int numValid;
  for (int ti = 0; ti < trackPts.size() ; ti++ ){ //index track
    cout<<"track="<<ti<<endl;
    InterpolateInvalidPoint(trackPts[ti]); 
    FindValidPoints(trackPts[ti], halfWindow, numValid);
    if (numValid > 0){
      ComputeGradient(trackPts[ti], halfWindow);
      ComputeSalientFeatures(trackPts[ti],  topPercent, numValid);
      MakeLinearWeights(trackPts[ti], halfWindow);

      stringstream ss;
      ss<<ti;
      string lolaTrackFeaturesFilename = prefix_from_filename(filename) + ss.str() +".txt";  
      SaveWeights(trackPts[ti], lolaTrackFeaturesFilename);
    } 
  }
}

//this function loops through each set of track points separatly, 
//if we see the following pattern of points (valid, in-valid, valid)
//the non-valid point will have it's reflectance set to the average of the two valid points
//
//above non-valid was defined as reflectance = -1 
int InterpolateInvalidPoint(vector< LOLAShot>& trackPts)
{
   for (int si = 1; si < trackPts.size()-1; si++ ){
   
      if(trackPts[si].reflectance < 0){
        if ((trackPts[si-1].reflectance > 0 ) && (trackPts[si+1].reflectance > 0)){
           trackPts[si].reflectance = (trackPts[si-1].reflectance + trackPts[si+1].reflectance )/2;  
        }
      }

    }
 
    return 0;
}



//function loops through points setting the LOLAShot integer field .calc_acp to "1" 
//if the point is valid & there are sufficent valid points on either side of that point compute a gradient.  
//otherwise .calc_acp is set to 0 if either the point is invalid or there are insufficient 
//valid points bordering that point to compute a gradient
//halfWindow = (length of the derivative filter-1)/2
//numValid = the number of points for which the derivative filter is computed
 
int FindValidPoints(vector<LOLAShot > & trackPts,int halfWindow, int &numValid)
{

  int t_inv = -1; //trailing invalid
  int l_inv = -1; //leading invalid
  numValid = 0;
 
  int si; //shot  index
  
  int numShots = trackPts.size();
  t_inv = -1;
  l_inv = -1;
  
  for (si = 0; si < numShots; si++){

      //don't compute the derivatives for points close to the track boundary
      if (si >= numShots-halfWindow){
         trackPts[si].calc_acp = false;
         continue;
      }

      if (trackPts[si].reflectance < 0){  
          t_inv = si;
      }

      if (trackPts[si+halfWindow].reflectance < 0){
          l_inv = si + halfWindow;
      }

      if ((si > (halfWindow-1)) && (si < (numShots-halfWindow+1)) ){
        
        if ( ((si-t_inv) > halfWindow) && (l_inv < si)){
           trackPts[si].calc_acp = 1; 
           numValid ++;
        }
        else{
           trackPts[si].calc_acp = 0;
        }

      }

      else{
          trackPts[si].calc_acp = 0;
      } 
    
  }
 
  cout<<"numValid in FindValidPts= "<<numValid<<endl;
  return 0;
}



//calculates the gradient of reflectance using a filter of shape: [ -1x--edge-- 0  +1x--edge--]
int ComputeGradient(vector<LOLAShot >& trackPts, int halfWindow)
{
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
  int numShots = trackPts.size();

  printf("trackPts.size() = %d\n",numShots);

  for(int si = halfWindow; si < (numShots - halfWindow); si++){

      if( trackPts[si].calc_acp){
        sum = 0;
        for(int place = 0; place < windowSize; place ++){
          sum += f[place]*trackPts[si-halfWindow+place].reflectance;
        }
        trackPts[si].filter_response = sum;
         //printf("track_pts[%d][%d] = %0.1f\n", iTrk , ip, track_pts[iTrk][ip].filter_response);
      }
  }
 
  return 0;
}


int ComputeSalientFeatures( vector< LOLAShot > & trackPts, float topPercent, int numValid)
{
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
  
  cout<<"numValid= "<<numValid<<endl;
  list_all_response.resize(numValid);

  for( int si = 0; si < trackPts.size(); si++ ){
 
      if(trackPts[si].calc_acp){
        list_all_response[numRead] = abs(trackPts[si].filter_response);
        numRead ++;
      }
  }
  
 
  // 1. sort list
  sort(list_all_response.begin(),list_all_response.end() );
 
 
  // 2. define the threshold
  int take_point = 0;
  take_point = (int)ceil( (1-topPercent) * numValid);
  //take_thresh = 0.0;
  float take_thresh = list_all_response[take_point];
  printf("list_all_response.size() = %d\n",(int)list_all_response.size());
  printf("make_weights report: take_point = %d, take_thresh = %f\n",take_point, take_thresh);
 

  // 3. assign values 
  int count_greater = 0;
  int others = 0;
 
  for(int si = 0; si < trackPts.size(); si ++ ){
    
      if( abs(trackPts[si].filter_response) >= take_thresh ){
        trackPts[si].weight_lsq = 1.0;  
        count_greater ++;
      }
      else{
        trackPts[si].weight_lsq = 0.05;
        others ++;
      }
  }

  printf("count_greater = %d, others = %d", count_greater, others);
  printf("Make weights function should be complete - take_point = %d, take_thresh = %f\n", take_point, take_thresh);
  return 0;
}




//assign linear weights (shapped like a pyramid) to all points around one of the interest point selected by make_weights
int MakeLinearWeights( vector< LOLAShot > & trackPts, const int &halfWindow)
{

  float windowSize_f = 2*halfWindow + 1;
  float halfWindow_f = halfWindow;

  float to_add = 0.0;
  printf("windowSize_f = %f, halfWindow_f = %f\n", windowSize_f, halfWindow_f);

    //zero the track before begining 
    for(int si = 0; si < trackPts.size(); si++){
      trackPts[si].weight_prd = 0; // just to begin
    }

    for( int si = 0; si < trackPts.size() ; si++){
      if( trackPts[si].weight_lsq == 1 ){

        //write out using edge - calculate normalization later
        to_add = 0.0;
        for( int m = 0; m < windowSize_f; m++){
          to_add = (windowSize_f - abs(m-halfWindow_f)) / windowSize_f;
          trackPts[si-halfWindow+m].weight_prd += to_add;
        }
      }
    }
 

  printf("Done computing the linear weights\n");
  return 0;
}


//write out all weight specific variables to a txt file specified by "f_name".  
//writen out is: reflectance, calc_acp, filter_response, weight_lsq, weight_prd
int SaveWeights(vector< LOLAShot>& trackPts, string filename)

{  
  cout << "saving to file "<<filename<<" ...";
 FILE* sFile;
 sFile = fopen(filename.c_str(), "w");

   for (int si = 0; si < trackPts.size(); si++ ){
     fprintf( sFile, "si= %d reflectance= %f calc_acp= %d filter_response= %f weight= %f weight_prd= %f \n", 
                     si, trackPts[si].reflectance, trackPts[si].calc_acp, 
                     trackPts[si].filter_response, trackPts[si].weight_lsq, trackPts[si].weight_prd );
   }

 fclose(sFile);
 cout << "done"<< endl;
}


