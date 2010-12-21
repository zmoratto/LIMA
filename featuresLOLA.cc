#include "featuresLOLA.h"


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
         trackPts[si].calc_acp = 0;
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
 
  cout<<"number of valid points for dreivative computation= "<<numValid<<endl;
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
      }
  }
 
  return 0;
}


int ComputeSalientFeatures( vector< LOLAShot > & trackPts, float topPercent, int numValid)
{
  // builds a list of the abs(gradient response), 


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
  
 
  //sort the list of LOLA points
  sort(list_all_response.begin(),list_all_response.end() );
 
 
  //determine the salient feature threshold
  int take_point = 0;
  take_point = (int)ceil( (1-topPercent) * numValid);
  float salientFeatureThresh = list_all_response[take_point];
  printf("ComputeSalientFeatures: salientFeatureThresh = %f\n", salientFeatureThresh);
 

  //determine the salient features 
  int numSalientFeatures = 0;
  int others = 0;
 
  for(int si = 0; si < trackPts.size(); si ++ ){
    
      if( abs(trackPts[si].filter_response) >= salientFeatureThresh ){
        trackPts[si].featurePtLOLA = 1;  
        numSalientFeatures ++;
      }
      else{
        trackPts[si].featurePtLOLA = 0;
        others ++;
      }
  }

  printf("numSalientFeatures = %d, others = %d", numSalientFeatures, others);

  return 0;
}





