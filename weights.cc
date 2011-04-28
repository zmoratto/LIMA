#include "featuresLOLA.h"
#include "weights.h"


//computes and saves the weights of a LOLA track
void ComputeWeights( vector< vector<LOLAShot> >& trackPts, int halfWindow, float topPercent, string filename ){
  
  int numValid;
  for (int ti = 0; ti < trackPts.size() ; ti++ ){ //index track
    //cout<<"track="<<ti<<endl;
    InterpolateInvalidPoint(trackPts[ti]); 
    FindValidPoints(trackPts[ti], halfWindow, numValid);
    cout<<"track="<<ti<<"numValid points = "<<numValid<<endl;
    if (numValid > 0){
      ComputeGradient(trackPts[ti], halfWindow);
      ComputeSalientFeatures(trackPts[ti],  topPercent, numValid);
      MakeLinearWeights(trackPts[ti], halfWindow);

      //stringstream ss;
      //ss<<ti;
      //string lolaTrackFeaturesFilename = prefix_from_filename(filename) + ss.str() +".txt";  
      //SaveWeights(trackPts[ti], lolaTrackFeaturesFilename);
    } 
  }
  
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
      trackPts[si].weightRefl = 0; // just to begin
    }

    for( int si = 0; si < trackPts.size() ; si++){
      if( trackPts[si].featurePtRefl == 1 ){

        //write out using edge - calculate normalization later
        to_add = 0.0;
        for( int m = 0; m < windowSize_f; m++){
          to_add = (windowSize_f - abs(m-halfWindow_f)) / windowSize_f;
          trackPts[si-halfWindow+m].weightRefl += to_add;
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
     fprintf( sFile, "si= %d reflectance= %f calc_acp= %d filter_response= %f weight_lsq= %f weightLOLA= %f \n", 
                     si, trackPts[si].reflectance, trackPts[si].calc_acp, 
                     trackPts[si].filter_response, trackPts[si].featurePtRefl, trackPts[si].weightRefl );
   }

 fclose(sFile);
 cout << "done"<< endl;
}


