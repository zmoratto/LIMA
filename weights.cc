#include "weights.h"

//this function loops through each set of track points seperatly, 
//if we see the following pattern of points (valid, in-valid, valid)
//the non-valid point will have it's reflectance set to the average of the two valid points
//
//above non-valid was defined as reflectance = -1 
int InterpolateInvalidPoint(vector<vector< LOLAShot> > & track_pts){

  for(int iTrk = 0;iTrk < track_pts.size() ; iTrk ++ ){ //index track

    //debugg inner loop - correct invalid points
    for(int iPt = 1; iPt < track_pts[iTrk].size()-1;iPt ++ ){
    
      //printf("reflectance at %d = %f\n",iPt,track_pts[iTrk][iPt].reflectance);
      if(track_pts[iTrk][iPt].reflectance < 0){
        if((track_pts[iTrk][iPt-1].reflectance > 0 ) && (track_pts[iTrk][iPt+1].reflectance > 0)){
          track_pts[iTrk][iPt].reflectance = (track_pts[iTrk][iPt-1].reflectance + track_pts[iTrk][iPt+1].reflectance )/2;  
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
int FindValidPoints(vector< vector<LOLAShot > > & track_pts,int edge, int & num_valid){


  int t_inv = -1;  //trailing invalid
  int l_inv = -1; //leading invalid
  num_valid = 0;
  for(int iTrk = 0; iTrk < track_pts.size(); iTrk ++){
    int n_shts = track_pts[iTrk].size();
    t_inv = -1;
    l_inv = -1;
    for(int ip = 0; ip < n_shts; ip++){
      //cout << "check_acceptable_computes: ip = " << ip;
      if(ip >= n_shts-edge){
        track_pts[iTrk][ip].calc_acp = false;
        //cout << "ip = " << ip << ", is false" << endl; 
        continue;
      }

      if((iTrk == 0)&&(ip<5)){
        printf("what the hell is going on?\n");
        printf("t_inv = %d, l_inv = %d, track_pts[%d][%d].reflectance = %f\n", t_inv, l_inv, iTrk, ip, track_pts[iTrk][ip].reflectance );//print out stuff see what's wrong 
      }

      if(track_pts[iTrk][ip].reflectance < 0){  
        t_inv = ip;
      }
      if(track_pts[iTrk][ip+edge].reflectance < 0){
        l_inv = ip + edge;
      }

      if((ip > (edge-1)) && (ip < (n_shts-edge+1)) ){
        if(ip==2){printf("Did we get here with the 2nd iteration?\n");}
        
        if( ((ip-t_inv) > edge) && (l_inv < ip)){
          track_pts[iTrk][ip].calc_acp = 1;
          num_valid ++;
          //cout <<"ip = " << ip << ", is true"  << endl;
        }else{
          track_pts[iTrk][ip].calc_acp = 0;
          //cout <<"ip = " << ip << ", is false"  << endl;
        }
      }else{
        track_pts[iTrk][ip].calc_acp = 0;
        //cout << "ip = " << ip << ", is false" << endl;
      } 
      /*
         if(track_pts[iTrk][ip].calc_acp ){
         printf("in function check: track_pts[0][%d].calc_acp = true\n", ip );
         } else{  
         printf("in function check: track_pts[0][%d].calc_acp = flase\n", ip );
         }
       */
    }
  }
  printf("\n\nEnd calc_accptable\n\n"); 
  return 0;
}

//2.int filter
//calculates the gradient of reflectance using a filter of shape: [ -1x--edge-- 0  +1x--edge--]
int ComputeGradient(vector< vector<LOLAShot > > & track_pts, int edge){
  
  cout << "Starting grad filter application\n" << endl;
  //build filter
  vector<float> f;
  float sum = 0; // filter response
  int width = 2*edge + 1;
  f.resize(width);
  for(int i = 0; i < width ;i++ ){
    if(i < edge){f[i] = -1.0; }
    else if(i==edge){f[i] = 0.0;}
    else{f[i] = 1.0;}
    //view the filter
    //printf("f[%d] = %f",i,f[i]);
  }
  cout << endl;

  //apply filter
  int iTrk = 0; // eventuall loop over all of this!
  //cout << "Get to start of filter application?" << endl;;
  for(int iTrk = 0; iTrk < track_pts.size(); iTrk++){
    int n_shts = track_pts[iTrk].size();
    printf("track_pts[%d].size() = %d\n",iTrk,n_shts);
    for(int ip = edge; ip < (n_shts - edge); ip++){
      if( track_pts[iTrk][ip].calc_acp){
        sum = 0;
        for(int place = 0; place < width; place ++){
          sum += f[place]*track_pts[iTrk][ip-edge+place].reflectance;
        }
        track_pts[iTrk][ip].filter_response = sum;
      
         //printf("track_pts[%d][%d] = %0.1f\n", iTrk , ip, track_pts[iTrk][ip].filter_response);
      }else{
       // printf("n. calc not acceptable, iTrk = %d, ip = %d\n", iTrk, ip);
      }
      // printf("track_pts[%d][%d] = %0.1f\n", iTrk , ip, track_pts[iTrk][ip].filter_response);
    }
  }
  cout << "End of filter application" << endl;
  return 0;
}


// we'll start with [0,1] based on top 5% 
// then allow room for growth
//3.int make_weights
int ComputeSalientFeatures( vector<vector< LOLAShot > > & track_pts, int edge, float take_p, int& num_valid, float& take_thresh){
  // builds a list of the abs(gradient response), 
  // finds the response value that corresponds with the top "take_p" %
  // assignes LOLAShot.weight_lsq = 1 if the response of that point is in the top take_p of all points
  // 
  // 
  //make weights - start with top 'take_p' accross all tracks
  //initial implementation will be slow

  // 0. make a list
  vector<float> list_all_response;
  int num_read = 0;
  //int iTrk = 0;

  list_all_response.resize(num_valid);
  for(int iTrk = 0; iTrk < track_pts.size(); iTrk++ ){
    //printf("list_all_response.size() = %d\n", list_all_response.size() );

    for( int i = 0; i < track_pts[iTrk].size(); i++ ){
      /*
         if(track_pts[iTrk][i].calc_acp ){
         printf("in function check: track_pts[0][%d].calc_acp = true\n", i );
         } else{  
         printf("in function check: track_pts[0][%d].calc_acp = flase\n", i );
         }*/
      // printf("track_pts[%d][%d].calc_acp = %f\n", iTrk, i, track_pts[iTrk][i].calc_acp );
      if(track_pts[iTrk][i].calc_acp){
        list_all_response[num_read] = abs(track_pts[iTrk][i].filter_response);
        num_read ++;
        //write out all entries after each input
        /*printf("we have just read the %dth entery into list_all_response:",num_read );
          for(int i_exam=0; i_exam<list_all_response.size(); i_exam++){
          printf(" %f\n", list_all_response[i_exam] );    
          }
         
        cout << endl;
      */
      }
    }
  } 
  //printf("num_valid = %d, num_read = %d\n ", num_read, num_valid );

  // 1. sort list
  sort(list_all_response.begin(),list_all_response.end() );
  //debugg by viewing this stuff
  //vector<float>::iterator it;  
  //for(it=list_all_response.begin(); it!=list_all_response.end(); it++){
  //  cout << " " << *it;}
  //cout << endl;
  /*
     printf("Try examining list_all_response again:\n");
     for(int t_read = 0; t_read < list_all_response.size(); t_read ++){
      printf("list_all_response[%d] = %f\n ", t_read, list_all_response[t_read]);
     }
   
  cout << endl;*/
  // 2. define cut
  int take_point = 0;
  take_point = (int)ceil( (1-take_p) * num_valid);
  take_thresh = 0.0;
  take_thresh = list_all_response[take_point];
  printf("list_all_response.size() = %d\n",list_all_response.size());
  printf("make_weights report: take_point = %d, take_thresh = %f\n",take_point, take_thresh);
  // 3. assign values - 
  //cout << "Do we get to the assignment of weights?" << endl;
  int count_greater = 0;
  int others = 0;
  for(int t_weight = 0; t_weight < track_pts.size();t_weight++ ){
    for(int itr_all = 0; itr_all < track_pts[t_weight].size(); itr_all ++ ){
      //cout << "itr_all = " << itr_all << endl;
      if( abs(track_pts[t_weight][itr_all].filter_response) >= take_thresh ){
        track_pts[t_weight][itr_all].weight_lsq = 1.0;  
        count_greater ++;
      }else{
        track_pts[t_weight][itr_all].weight_lsq = 0.05;
        others ++;
      }
      //cout << "Completed a loop - assign values" << endl;
      //printf("itr_all = %d, track_pts[%d][%d]  = %f\n", itr_all, t_weight, itr_all, track_pts[t_weight][itr_all].weight_lsq );
    }
  }
  printf("count_greater = %d, others = %d", count_greater, others);
  printf("Make weights function should be complete - take_point = %d, take_thresh = %f\n", take_point, take_thresh);
  return 0;
}

//assign linear weights (shapped like a pyramid) to all points around one of the interest point selected by make_weights
int MakeLinearWeights( vector< vector< LOLAShot > > & track_pts, const int &edge, const int & width){

  float w_f = width;
  float e_f = edge;
  float to_add = 0.0;
  printf("w_f = %f, e_f = %f\n", w_f, e_f);
  for(int iTrk = 0; iTrk < track_pts.size(); iTrk++ ){
    //zero the track before begining 
    for(int IP = 0; IP < track_pts[iTrk].size(); IP++){
      track_pts[iTrk][IP].weight_prd = 0; // just to begin
    }

    for( int ip = 0; ip < track_pts[iTrk].size() ; ip++){
      if( track_pts[iTrk][ip].weight_lsq == 1 ){
        //write out using edge - calculate normalization latter
        to_add = 0.0;
        for( int m = 0; m < width; m++){
          to_add = (w_f - abs(m-e_f)) / w_f;
          track_pts[iTrk][ip-edge+m].weight_prd += to_add;
         // if(iTrk >= 2 ){
         //    printf("m = %d, ip = %d, (ip-edge+m) = %d, track_pts[%d][%d].weight_prd = %f, to_add = %f \n", m, ip, (ip-edge+m ), iTrk, (ip-edge+m), track_pts[iTrk][ip].weight_prd, to_add);
         //  }
        }
      }
    }
  }
  printf("Done computing the linear weights\n");
  return 0;
}

//write out all weight specific variables to a txt file specified by "f_name".  
//writen out is: reflectance, calc_acp, filter_response, weight_lsq, weight_prd
int SaveWeights( vector< vector< LOLAShot> > & track_pts, string & f_name){
  
 cout << "Begining the end - write out commensing " << endl;
 FILE* sFile;
 sFile = fopen(f_name.c_str(), "w");
 fprintf(sFile, "Write weights:\n");
 for (int iTrk = 0; iTrk < track_pts.size(); iTrk++ ){
   for (int ip = 0; ip < track_pts[iTrk].size(); ip++ ){
     fprintf( sFile, "iTrk=  %d ip= %d reflectance= %f calc_acp= %d filter_response= %f weight= %f weight_prd= %f \n", 
                     iTrk, ip, track_pts[iTrk][ip].reflectance, track_pts[iTrk][ip].calc_acp, 
                     track_pts[iTrk][ip].filter_response, track_pts[iTrk][ip].weight_lsq, track_pts[iTrk][ip].weight_prd );
   }
 }

 fclose(sFile);
 cout << "Write out completed"<< endl;
}

/*
int main(){
  // debugging! - arrays with test cases
  vector<float> refl_base;
  refl_base.resize(14);
  refl_base[0] = -1;
  refl_base[1] = 2;
  refl_base[2] = 3;
  refl_base[3] = -1; 
  refl_base[4] = 5;
  refl_base[5] = 6;
  refl_base[6] = 7;
  refl_base[7] = 8;
  refl_base[8] = -1;
  refl_base[9] = 10;
  refl_base[10] = 11;
  refl_base[11] = 14;
  refl_base[12] = -1;
  refl_base[13] = -1;

 */ /*
     refl_base.resize(14);
     refl_base[0] = -1;
     refl_base[1] = 2;
     refl_base[2] = 9;
     refl_base[3] = -1; 
     refl_base[4] = 5;
     refl_base[5] = 6;
     refl_base[6] = 7;
     refl_base[7] = 8;
     refl_base[8] = -1;
     refl_base[9] = 1;
     refl_base[10] = 11;
     refl_base[11] = 12;
     refl_base[12] = -1;
   */
/*
  vector<float> refl_two;
  refl_two.resize(14);
  refl_two[0] = 5;
  refl_two[1] = -1;
  refl_two[2] = 8;
  refl_two[3] = 11;
  refl_two[4] = 10;
  refl_two[5] = 12;
  refl_two[6] = 11;
  refl_two[7] = -1;
  refl_two[8] = 15;
  refl_two[9] = 9;
  refl_two[10] = 20;
  refl_two[11] = 10;
  refl_two[12] = -1;
  refl_two[13] = 17;


  vector<float> refl_3;
  refl_3.resize(14);
  refl_3[0] = -1;
  refl_3[1] = -1;
  refl_3[2] = -1;
  refl_3[3] = -1;
  refl_3[4] = 0.1;
  refl_3[5] = 0.5;
  refl_3[6] =  1;
  refl_3[7] = 5;
  refl_3[8] = 10;
  refl_3[9] = -1;
  refl_3[10] = -1;
  refl_3[11] = -1;
  refl_3[12] = -1;
  refl_3[13] = -1;


  vector<float> refl_4;
  refl_4.resize(14);
  refl_4[0] = -1;
  refl_4[1] = 1;
  refl_4[2] = 101;
  refl_4[3] = 102;
  refl_4[4] = 99;
  refl_4[5] = 99;
  refl_4[6] = 102;
  refl_4[7] = 101;
  refl_4[8] = 1;
  refl_4[9] = -1;
  refl_4[10] = -1;
  refl_4[11] = -1;
  refl_4[12] = -1;
  refl_4[13] = -1;


  // read in reflectance values
  vector< vector< LOLAShot > > track_pts; 
  track_pts.resize(4);
  track_pts[0].resize(14);
  track_pts[1].resize(14);
  track_pts[2].resize(14);
  track_pts[3].resize(14);
  int num_tracks = track_pts.size();
  cout << "Main function initialized reflectance: "<< endl;
  for(int k = 0; k < 4; k++){
    int i = 0;
    for(i; i < refl_base.size(); i++){
      if(k ==0 ){
        track_pts[k][i].reflectance = refl_base[i];
        cout << "track_pts[k][" << i <<"]= " << track_pts[k][i].reflectance << endl;;
      }else if(k == 1){
        cout << endl;
        track_pts[k][i].reflectance = refl_two[i];
        printf("track_pts[%d][%d].reflectance = %f", k, i, track_pts[k][i].reflectance);
      }else if(k==2){
        cout << endl;
        track_pts[k][i].reflectance = refl_3[i];
        printf("track_pts[%d][%d].reflectance = %f", k, i, track_pts[k][i].reflectance);
      }else if(k==3){
        cout << endl;
        track_pts[k][i].reflectance = refl_4[i];
        printf("track_pts[%d][%d].reflectance = %f", k, i, track_pts[k][i].reflectance);
      }
    }
  }

  // linear interpolation of invalid points
  correct_invalid_linear( track_pts);

  cout << "Main function initialized reflectance: "<< endl;
  for(int K = 0; K < num_tracks; K++ ){
    for(int i_after = 0; i_after < refl_base.size(); i_after++){
      cout << "track_pts["<< K << "][" << i_after <<"]= " << track_pts[K][i_after].reflectance << endl;;
    }
  }
  // check accepatble computes
  int edge = 2;
  int num_valid = 0;
  check_acceptable_computes( track_pts, edge, num_valid);
  printf("Main function acceptable computes check is: \n");
  int ITAR = 0;
  for(ITAR; ITAR < num_tracks; ITAR ++){
    for(int i_p_filter = 0; i_p_filter < refl_base.size(); i_p_filter++){
      if(track_pts[ITAR][i_p_filter].calc_acp){
        printf("track_pts[%d][%d].calc_acp = true\n", ITAR, i_p_filter );
      } else{
        printf("track_pts[%d][%d].calc_acp = false\n", ITAR, i_p_filter );
      }
    } 
  }

  cout << endl;

  // filter
  grad_filter( track_pts, edge);
  float take_p = 0.2; // change from debugg to actual examples
  float take_thresh = 0.0;
  // decide weights
  cout << "Starting: make_weights_not_war " << endl;
  make_weights_not_war( track_pts, edge, take_p, num_valid, take_thresh);
  // view weights outside of the internal function
  printf("Print out weights_lsq after assignment: ");
  ITAR = 0;
  int i_prt = 0;
  for(ITAR; ITAR<num_tracks; ITAR ++){
    for(i_prt = 0; i_prt < track_pts[ITAR].size(); i_prt ++){
      printf("track_pts[%d][%d].weight_lsq =  %f\n",ITAR, i_prt, track_pts[ITAR][i_prt].weight_lsq);
    }
  }
  cout << endl << endl << "The End."  << endl;

  // compute test cases - check

  // write out the results
  string f_name = "../results/weights_record.txt";
  pyramid_weights( track_pts, edge, 2*edge+1);
  ITAR = 0;
  i_prt = 0;
  for( ITAR; ITAR<num_tracks; ITAR ++){
    for(i_prt = 0; i_prt < track_pts[ITAR].size();i_prt++ ){
      printf("track_pts[%d][%d].weight_prd =  %f\n",ITAR, i_prt, track_pts[ITAR][i_prt].weight_prd);
    }
  }
  
  write_refl_and_weights( track_pts, f_name);
  return 0;
}
*/
