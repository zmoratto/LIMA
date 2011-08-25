// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <boost/tokenizer.hpp>
#include <vw/Core.h>
#include <vw/Math/Matrix.h>
#include <vw/Cartography.h>
#include <vw/Cartography/SimplePointImageManipulation.h>
#include <vw/FileIO.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include "match.h"
#include "coregister.h"
#include "display.h"
#include "util.h"
#include "tracks.h"

using namespace vw;
using namespace vw::cartography;
using namespace vw::math;

//needed for map to cam back projection - START
#include <Camera.h>
#include <Pvl.h>
#include <CameraFactory.h>
#include <CameraDetectorMap.h>
//needed for map to cam back projection - END

//Constructor for pointCloud which is really just a Vector3 
// with some extra data fields.
pointCloud::pointCloud
  (
  Vector3 coords,
  int     y,
  int     mo,
  int     d,
  int     h,
  int     mi,
  float   s,
  int     detector 
  )
  :
  Vector<double,3>( coords )
  {
  year = y;
  month = mo;
  day = d;
  hour = h;
  min = mi;
  sec = s;
  s = detector;
  }

void
LOLAShot::init
  (
  vector<pointCloud> pc,
  vector<imgPoint>   ip,
  vector<DEMPoint>   dp,
  int                v,
  int                cpi,
  float              r,
  float              si,
  int                ca,
  float              fr,
  float              fpr,
  float              wr,
  float              fpL,
  float              fL
  )
  {
  LOLAPt = pc;
  imgPt = ip;
  DEMPt = dp;
  valid = v;
  centerPtIndex = cpi;
  reflectance = r;
  synthImage = si;
  calc_acp = ca;
  filter_response = fr;
  featurePtRefl = fpr;
  weightRefl = wr;
  featurePtLOLA = fpL;
  filresLOLA = fL;
  }

LOLAShot::LOLAShot( vector<pointCloud> pcv )
  {
  LOLAShot::init( pcv );
  }
LOLAShot::LOLAShot( pointCloud pt )
  {
  vector<pointCloud> pcv( 1, pt );
  LOLAShot::init( pcv );
  }



vector<vector<LOLAShot> > CSVFileRead(string CSVFilename)
	{
	// This is specifically for reading the RDR_*PointPerRow_csv_table.csv files only.
	// It will ignore lines which do not start with year (or a value that can be converted
	// into an integer greater than zero, specfically).
	ifstream myfile (CSVFilename.c_str());
 
	int trackIndex = 0;
	vector<vector<LOLAShot> > trackPts(1);
 
	pointCloud currPt;
	pointCloud prevPt;

	if (!myfile)
		{
		vw_throw( vw::IOErr() << "Unable to open track file \"" << CSVFilename << "\"" );
		}

	// Prepare the tokenizer stuff.
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(" ,");

	// myfile >> ignoreLine; // Skip header line

	string line;
    while ( getline( myfile, line, '\n' ) )
		{

		/*
		// For Debugging: this reads the first line, parses, numbers, and prints
		tokenizer tokens( line, sep );
		int counter(0);
		for( tokenizer::iterator token = tokens.begin();
				token != tokens.end();
				++token, counter++)
			{
			cerr	<< boost::lexical_cast<std::string>(counter) 
					<< ". " << *token << endl;
			}
		exit(1);
		*/

		tokenizer tokens( line, sep );
		tokenizer::iterator token = tokens.begin();

		istringstream date_s( *token );
		date_s >>  setw(4) >> currPt.year
				>> ignoreOne // -
				>> setw(2) >> currPt.month
				>> ignoreOne // -
				>> setw(2) >> currPt.day
				>> ignoreOne // T
				>> setw(2) >> currPt.hour
				>> ignoreOne // :
				>> setw(2) >> currPt.min
				>> ignoreOne // :
				>> currPt.sec;
		//cout<<"year "<< currPt.year << endl;
		if( currPt.year <= 0 ) { continue; }
		//cout << currPt.year <<' '<< currPt.month <<' '<< currPt.day 
		//	<<' '<< currPt.hour <<' '<< currPt.min <<' '<< currPt.sec << endl;

		++token;
		istringstream lon_s( *token );
		lon_s >> currPt.x(); // Pt_Longitude
		// cout<<"lon "<< currPt.coords	s_s >> "s_s"<<currPt.s; // S(0) << endl;

		++token;
		istringstream lat_s( *token );
		lat_s >> currPt.y(); // Pt_Latitude

		++token;
		istringstream rad_s( *token );
		rad_s >> currPt.z(); // Pt_Radius
		// cout<<"radius "<< currPt.coords(2) << endl;

		advance( token, 8 );
		istringstream s_s( *token );
		s_s >>currPt.s; // S
                //cout<<"s_s "<<currPt.s<<endl;

        if ((currPt.x()!=0.0) && (currPt.y()!=0.0) )
			{ //valid lidar point

			if( trackPts[trackIndex].empty() ) {
			  trackPts[trackIndex].push_back( LOLAShot(currPt) );
            }
			else
				{
				if( GetTimeDiff(prevPt, currPt, 3000) )
					{ //new track
					trackIndex++;
					trackPts.resize(trackIndex+1); //start new track
					trackPts[trackIndex].push_back( LOLAShot(currPt) ); //put the new shot in it
//add last shot to the previous track
					//shot.LOLAPt.clear();
					//trackIndex++;
					//trackPts.resize(trackIndex+1); //start new track
					//shot.LOLAPt.push_back(currPt);
	    		 	}
				else
					{ //same track
					if( GetTimeDiff(prevPt, currPt, 0) )
						{//new shot
						trackPts[trackIndex].push_back( LOLAShot(currPt) );
						//shot.LOLAPt.clear();
						//shot.LOLAPt.push_back(currPt);
						}
					else
						{ //same shot
						// shot.LOLAPt.push_back(currPt);
						trackPts[trackIndex].back().LOLAPt.push_back( currPt );
						}
					}
				}
    
			//copy current pc into prevPt
			prevPt = currPt;
            /*
			prevPt.x() = currPt.x();
			prevPt.y() = currPt.y();
			prevPt.z() = currPt.z();
			prevPt.year = currPt.year;
			prevPt.month = currPt.month;
			prevPt.day = currPt.day;
			prevPt.hour = currPt.hour;
			prevPt.min = currPt.min;
			prevPt.sec = currPt.sec;   
			prevPt.s = currPt.s;
			*/
			}

		} 

	myfile.close();

	return trackPts; 
	}

Vector2 ComputeMinMaxValuesFromCub(string cubFilename)
{
  Vector2 minmax;
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodataVal = rsrc->nodata_read();
  //cout<<"nodaval:"<<nodataVal<<endl;
  DiskImageView<PixelGray<float> > isis_view( rsrc );
  int width = isis_view.cols();
  int height = isis_view.rows();
  float minVal = 100000000.0;
  float maxVal = -100000000.0;
  for (int i = 0; i < height; i++){
    for (int j = 0; j < width; j++){
      
      if ((isis_view(j,i) < minVal) && (isis_view(j,i) > nodataVal)){
	minVal = isis_view(j,i);
      }
      if ((isis_view(j,i) > maxVal) && (isis_view(j,i) > nodataVal)){
	maxVal = isis_view(j,i);
      }
    }
  }
  minmax(0) = minVal;
  minmax(1) = maxVal;
  cout<<"min="<<minVal<<", max="<<maxVal<<endl;

  return minmax;
}

int GetAllPtsFromCub(vector<vector<LOLAShot > > &trackPts, string cubFilename)
{

  vector<pointCloud> LOLAPts;
  
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodata_value = rsrc->nodata_read();
  cout<<"no_data_val="<<nodata_value<<endl; 

  DiskImageView<PixelGray<float> > isis_view( rsrc );
  int width = isis_view.cols();
  int height = isis_view.rows();
  camera::IsisCameraModel model(cubFilename);

  //calculate the min max value of the image
  Vector2 minmax = ComputeMinMaxValuesFromCub(cubFilename);

  float minTrackVal =  1000000.0;
  float maxTrackVal = -1000000.0;

  int numValidImgPts = 0;

  InterpolationView<EdgeExtensionView<DiskImageView<PixelGray<float> >, ConstantEdgeExtension>, BilinearInterpolation> interpImg
    = interpolate(isis_view, BilinearInterpolation(), ConstantEdgeExtension());
  
  for(unsigned int k = 0; k < trackPts.size();k++){
    for(unsigned int i = 0; i < trackPts[k].size(); i++){
      
      trackPts[k][i].valid = 1; 

      LOLAPts = trackPts[k][i].LOLAPt;
  
      trackPts[k][i].imgPt.resize(LOLAPts.size());

      for (unsigned int j = 0; j < LOLAPts.size(); j++){
          
	    float lon = LOLAPts[j].x();
	    float lat = LOLAPts[j].y();
	    float rad = LOLAPts[j].z();
            
            Vector3 lon_lat_rad (lon,lat,rad*1000);
            Vector3 xyz = cartography::lon_lat_radius_to_xyz(lon_lat_rad);
            Vector2 cub_pix = model.point_to_pixel(xyz);
            float x = cub_pix[0];
            float y = cub_pix[1];
            //check that (x,y) are within the image boundaries
	    if ((x>=0) && (y>=0) && (x<width) && (y<height)){//valid position  
              //check for valid data as well
              if (interpImg(x, y)>minmax(0)){//valid values
		 trackPts[k][i].imgPt[j].val = interpImg(x, y);
	         trackPts[k][i].imgPt[j].x = cub_pix[0];
	         trackPts[k][i].imgPt[j].y = cub_pix[1];
                 if (interpImg(x,y) < minTrackVal){
		     minTrackVal = interpImg(x,y);
                 }
		 if (interpImg(x,y) > maxTrackVal){
		     maxTrackVal = interpImg(x,y);
                 }
                 numValidImgPts++;
	      }
              else{//invalidate the point
                 trackPts[k][i].valid = 0;
              }
	    }
            else{ //invalidate the point  
                 trackPts[k][i].valid = 0; 
            }
	
      }
      /*
      // this must go to computeReflectance function - START
      //check for valid shot with 3 points to compute reflectance
      pointCloud centerPt  = GetPointFromIndex( LOLAPts, 1);
      pointCloud topPt     = GetPointFromIndex( LOLAPts, 3);
      pointCloud leftPt    = GetPointFromIndex( LOLAPts, 2);
      if ((centerPt.s == -1) || (topPt.s == -1) || (leftPt.s == -1) || (LOLAPts.size() > 5)){//invalid LOLA shot
          trackPts[k][i].valid = 0; 
      }
      // this must go to computeReflectance function - END
      */
    }//i  
  }//k
  cout <<"minTrackVal="<<minTrackVal<<", maxTrackVal="<<maxTrackVal<<endl;
  
  return numValidImgPts;
}

vector<float> GetTrackPtsByID(vector<LOLAShot> trackPts, int ID)
{
  vector<float> pts;
  for (unsigned int i = 0; i < trackPts.size(); i++){
    for (unsigned int k = 0; k < trackPts[i].LOLAPt.size(); k++){

      float rad = trackPts[i].LOLAPt[k].z();
      int id = trackPts[i].LOLAPt[k].s;

      if (id == ID){
        pts.push_back(rad);
      }

    }
  }

  return pts;
}

//computes the scale factor for all tracks at once
float ComputeScaleFactor(vector<vector<LOLAShot > >&trackPts)
{
  float nominator = 0.0;
  int numValidPts = 0;
  float scaleFactor;// = 1;

  //Vector2 minmax = ComputeMinMaxValuesFromCub(string cubFilename);

  for(unsigned int k = 0; k < trackPts.size();k++){
    for(unsigned int i = 0; i < trackPts[k].size(); i++){
      if ((trackPts[k][i].valid == 1) && (trackPts[k][i].reflectance != 0) &&(trackPts[k][i].reflectance != -1)){//valid track and non-zero reflectance

        //update the nominator for the center point

        for (unsigned int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){
          if (trackPts[k][i].LOLAPt[j].s == 1){
            //cout<< trackPts[k][i].imgPt[j].val/trackPts[k][i].reflectance<<endl;
            nominator = nominator + (trackPts[k][i].imgPt[j].val/*-minmax(0)*/)/trackPts[k][i].reflectance;
            numValidPts++;
          }
        }

        //update the denominator
        //numValidPts++;
      }
    }
  }

  //cout<<"NUM_VALID_POINTS="<<numValidPts<<endl;
  if (numValidPts != 0){ 
    scaleFactor = nominator/numValidPts;
  }
  else{
    //invalid scaleFactor, all tracks are invalid
    scaleFactor = -1;
  }
  return scaleFactor;
}

//computes the scale factor for all tracks at once
Vector2 ComputeGainBiasFactor(vector<vector<LOLAShot > >&trackPts)
{
  int numValidPts = 0;
  float sum_rfl = 0.0; 
  float sum_img = 0.0;
  float sum_rfl_2 = 0.0;
  float sum_rfl_img = 0.0;
  Matrix<float,2,2> rhs;
  Vector<float,2> lhs;
  Vector2 gain_bias;

  for(unsigned int k = 0; k < trackPts.size();k++){
    for(unsigned int i = 0; i < trackPts[k].size(); i++){
      if ((trackPts[k][i].valid == 1) && (trackPts[k][i].reflectance != 0) &&(trackPts[k][i].reflectance != -1)){//valid track and non-zero reflectance

        //update the nominator for the center point

        for (unsigned int j = 0; j < trackPts[k][i].LOLAPt.size(); j++){
          if (trackPts[k][i].LOLAPt[j].s == 1){
            sum_rfl = sum_rfl + trackPts[k][i].reflectance;
            sum_rfl_2 = sum_rfl_2 + trackPts[k][i].reflectance*trackPts[k][i].reflectance;
            sum_rfl_img = sum_rfl_img + trackPts[k][i].reflectance*(trackPts[k][i].imgPt[j].val);
            sum_img = sum_img + (trackPts[k][i].imgPt[j].val);
            numValidPts++;
          }
        }

        //update the denominator
        //numValidPts++;
      }
    }
  }

  //cout<<"NUM_VALID_POINTS="<<numValidPts<<endl;
  if (numValidPts != 0){ 
    rhs(0,0) = sum_rfl_2;
    rhs(0,1) = sum_rfl;
    rhs(1,0) = sum_rfl;
    rhs(1,1) = numValidPts;
    lhs(0) = sum_rfl_img;
    lhs(1) = sum_img;
    solve_symmetric_nocopy(rhs,lhs);
    gain_bias = lhs;
    cout<<"gain and bias"<<lhs<<endl;
  }
  else{
    //invalid scaleFactor, all tracks are invalid
    gain_bias(0) = 0.0;
    gain_bias(1) = 0.0;
  }
  return gain_bias;
}

Vector3 ComputePlaneNormalFrom3DPoints(vector<Vector3> pointArray)
{
  Vector3 normal;
  Matrix<float,5,4> rhs;
  Vector<float,3> lhs;
  for (unsigned int i = 0; i < pointArray.size(); i++){
    rhs(i,0)=pointArray[i][0];
    rhs(i,1)=pointArray[i][1]; 
    rhs(i,2)=pointArray[i][2];
    rhs(i,3)=1;
  }
  
  solve_symmetric_nocopy(rhs,lhs);
  return normal;

}

Vector3 ComputeNormalFrom3DPointsGeneral(Vector3 p1, Vector3 p2, Vector3 p3) {
  return -normalize(cross_prod(p2-p1,p3-p1));
}

float ComputeLunarLambertianReflectanceFromNormal(Vector3 sunPos, Vector3 viewPos, Vector3 xyz, Vector3 normal) 
{
  float reflectance;
  float L;

  //compute /mu_0 = cosine of the angle between the light direction and the surface normal.
  //sun coordinates relative to the xyz point on the Moon surface
  //Vector3 sunDirection = -normalize(sunPos-xyz);
  Vector3 sunDirection = normalize(sunPos-xyz);
  float mu_0 = dot_prod(sunDirection,normal);

  //compute  /mu = cosine of the angle between the viewer direction and the surface normal.
  //viewer coordinates relative to the xyz point on the Moon surface
  Vector3 viewDirection = normalize(viewPos-xyz);
  float mu = dot_prod(viewDirection,normal);

  //compute the phase angle /alpha between the viewing direction and the light source direction
  float rad_alpha, deg_alpha;
  float cos_alpha;

  cos_alpha = dot_prod(sunDirection,viewDirection);
  if ((cos_alpha > 1)||(cos_alpha< -1)){
    printf("cos_alpha error\n");
  }

  rad_alpha = acos(cos_alpha);
  deg_alpha = rad_alpha*180.0/M_PI;

  //printf("deg_alpha = %f\n", deg_alpha);

  //Bob Gaskell's model
  //L = exp(-deg_alpha/60.0);

#if 0 // trey
  // perfectly valid for alpha to be greater than 90?
  if (deg_alpha > 90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
  if (deg_alpha < -90){
    //printf("Error!!: rad_alpha = %f, deg_alpha = %f\n", rad_alpha, deg_alpha);
    return(0.0);
  }
#endif

  //Alfred McEwen's model
  float A = -0.019;
  float B =  0.000242;//0.242*1e-3;
  float C = -0.00000146;//-1.46*1e-6;

  L = 1.0 + A*deg_alpha + B*deg_alpha*deg_alpha + C*deg_alpha*deg_alpha*deg_alpha;

  //        std::cout << " sun direction " << sunDirection << " view direction " << viewDirection << " normal " << normal;
  //        std::cout << " cos_alpha " << cos_alpha << " incident " << mu_0 << " emission " << mu;
  //printf(" deg_alpha = %f, L = %f\n", deg_alpha, L);

  //if (mu_0 < 0.15){ //incidence angle is close to 90 deg
  if (mu_0 < 0.0){
    //mu_0 = 0.15;
    return (0.0);
  }

  if (mu < 0.0){ //emission angle is > 90
    mu = 0.0;
    //return (0.0);
  }

  if (mu_0 + mu == 0){
    //printf("negative reflectance\n");
    reflectance = 0.0;
  }
  else{
    reflectance = 2*L*mu_0/(mu_0+mu) + (1-L)*mu_0;
  }
  if (reflectance < 0){
    //printf("negative reflectance\n");
    reflectance = 0;
  }
  return reflectance;
}

int  ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks,  Vector3 cameraPosition, Vector3 lightPosition, CoregistrationParams coregistrationParams)
{

  vector<pointCloud> LOLAPts;
  int numValidReflPts = 0;

  float minReflectance =  10000.0;
  float maxReflectance = -10000.0;
 
  for (unsigned int k = 0; k < allTracks.size(); k++ ){
    for (unsigned int i = 0; i < allTracks[k].size(); i++){
      LOLAPts = allTracks[k][i].LOLAPt;
      
      /*
      //old
      pointCloud centerPt = GetPointFromIndex(LOLAPts, 3);
      pointCloud topPt = GetPointFromIndex(LOLAPts, 2);
      pointCloud leftPt = GetPointFromIndex(LOLAPts, 1);
      */
      
      //new
      pointCloud centerPt = GetPointFromIndex(LOLAPts, 1);
      pointCloud topPt = GetPointFromIndex(LOLAPts, 3);
      pointCloud leftPt = GetPointFromIndex(LOLAPts, 2);
      
      if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1) /*&& (allTracks[k][i].valid == 1)*/ && ((LOLAPts.size() <= 5))){

        centerPt.z() = centerPt.z()*1000;
        topPt.z() = topPt.z()*1000;
        leftPt.z() = leftPt.z()*1000;
        
        Vector3 xyz = lon_lat_radius_to_xyz(centerPt);
        Vector3 xyzTop = lon_lat_radius_to_xyz(topPt);
        Vector3 xyzLeft = lon_lat_radius_to_xyz(leftPt);

   
        Vector3 normal = ComputeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);

        allTracks[k][i].reflectance = ComputeLunarLambertianReflectanceFromNormal(lightPosition, cameraPosition, xyz, normal); 
        
        if (allTracks[k][i].reflectance != 0){ 
	  if (allTracks[k][i].reflectance < minReflectance){
	    minReflectance = allTracks[k][i].reflectance;
	  } 
	  if (allTracks[k][i].reflectance > maxReflectance){
	    maxReflectance = allTracks[k][i].reflectance;
	  }
          numValidReflPts++;
	}        
      }
      else{
        allTracks[k][i].reflectance = -1;
      }
    }//i
  }//k
 
  cout<<"NEW:"<<"minReflectance="<<minReflectance<<", maxReflectance="<<maxReflectance<<endl;

  return numValidReflPts;
}

int  ComputeAllReflectance( vector< vector<LOLAShot> >  &allTracks, ModelParams modelParams,  CoregistrationParams coregistrationParams)
{

  vector<pointCloud> LOLAPts;
  int numValidReflPts = 0;

  GlobalParams globalParams;
  globalParams.reflectanceType = coregistrationParams.reflectanceType;
  globalParams.slopeType = 1;
  globalParams.shadowThresh = 40;
  globalParams.albedoInitType = 1;
  globalParams.exposureInitType = 1;

  float minReflectance =  10000.0;
  float maxReflectance = -10000.0;
 
  for (unsigned int k = 0; k < allTracks.size(); k++ ){
    for (unsigned int i = 0; i < allTracks[k].size(); i++){
      LOLAPts = allTracks[k][i].LOLAPt;
      /*
      pointCloud centerPt = GetPointFromIndex(LOLAPts, 3);
      pointCloud topPt = GetPointFromIndex(LOLAPts, 2);
      pointCloud leftPt = GetPointFromIndex(LOLAPts, 1);
      */
      pointCloud centerPt = GetPointFromIndex(LOLAPts, 1);
      pointCloud topPt = GetPointFromIndex(LOLAPts, 3);
      pointCloud leftPt = GetPointFromIndex(LOLAPts, 2);

      if ((centerPt.s != -1) && (topPt.s != -1) && (leftPt.s != -1) /*&& (allTracks[k][i].valid == 1)*/ && ((LOLAPts.size() <= 5))){

        Datum moon;
        moon.set_well_known_datum("D_MOON");


        centerPt.z() = (centerPt.z()-1737.4)*1000;
        topPt.z() = (topPt.z()-1737.4)*1000;
        leftPt.z() = (leftPt.z()-1737.4)*1000;

        Vector3 xyz = moon.geodetic_to_cartesian(centerPt);
        Vector3 xyzTop = moon.geodetic_to_cartesian(topPt);
        Vector3 xyzLeft = moon.geodetic_to_cartesian(leftPt);
        Vector3 normal = computeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);

        allTracks[k][i].reflectance = ComputeReflectance(normal, xyz, modelParams, globalParams);
        if (allTracks[k][i].reflectance != 0){ 
	  if (allTracks[k][i].reflectance < minReflectance){
	    minReflectance = allTracks[k][i].reflectance;
	  } 
	  if (allTracks[k][i].reflectance > maxReflectance){
	    maxReflectance = allTracks[k][i].reflectance;
	  }
          numValidReflPts++;
	}
        
      }
      else{
        allTracks[k][i].reflectance = -1;
      }
    }//i
  }//k
 
  cout<<"minReflectance="<<minReflectance<<", maxReflectance="<<maxReflectance<<endl;

  return numValidReflPts;
}

pointCloud GetPointFromIndex(vector<pointCloud> const &  LOLAPts, int index)
{
  pointCloud pt;
  pt.s = -1;//invalid pointCloud
  for(unsigned int i = 0;i < LOLAPts.size(); i++){
    if (LOLAPts[i].s == index){
      return LOLAPts[i];
    }
  }

  return pt;
}


//this function will be similar to GetAllPtsFromImage: 
//TO DO: rename to GetAllPtsFromDEM(vector<vector<LOLAShot > > &trackPts,  ImageViewBase<ViewT> const& DEM, GeoReference const &DEMGeo)
vector<float> GetTrackPtsFromDEM(vector<LOLAShot> trackPts, string DEMFilename, int ID)
{
  DiskImageView<PixelGray<float> >   DEM(DEMFilename);
  GeoReference DEMGeo;
  read_georeference(DEMGeo, DEMFilename);

  vector<float> demPts;
  demPts.resize(trackPts.size());

  ImageViewRef<PixelGray<float>  >  interpDEM = interpolate(edge_extend(DEM.impl(),
        ConstantEdgeExtension()),
      BilinearInterpolation());

  int index = 0;
  for (unsigned int i = 0; i < trackPts.size(); i++){
    demPts[index] = -1;
    for(unsigned int j = 0; j < trackPts[i].LOLAPt.size(); j++){

      float lon = trackPts[i].LOLAPt[j].x();
      float lat = trackPts[i].LOLAPt[j].y();
      //float rad = trackPts[i].LOLAPt[j].coords[2];
      int id = trackPts[i].LOLAPt[j].s;

      Vector2 DEM_lonlat(lon, lat);
      Vector2 DEM_pix = DEMGeo.lonlat_to_pixel(DEM_lonlat);

      int x = (int)DEM_pix[0];
      int y = (int)DEM_pix[1];

      PixelGray<float>  DEMVal = interpDEM(x, y);

      if ((id == ID) && (trackPts[i].valid)){
        //demPts.push_back((float)DEMVal);
        demPts[index] = (float)DEMVal; 
      }             
    }
    index++;
  }

  return demPts;
}

void SaveDEMPoints(vector< vector<LOLAShot> > &trackPts, string DEMFilename, string filename)

{    
  for (unsigned int k = 0; k < trackPts.size(); k++){
    vector<float> demPts = GetTrackPtsFromDEM(trackPts[k], DEMFilename, 3);
    string prefixTrackFilename =  prefix_from_filename(filename);  
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    SaveVectorToFile(demPts, string(trackFilename));     
  }
}


void SaveReflectancePoints(vector< vector<LOLAShot> >  &trackPts, Vector2 gain_bias, string filename)
{

  FILE *fp;
  float gain = gain_bias(0);
  float bias = gain_bias(1);
 
  for (unsigned int k = 0; k < trackPts.size(); k++ ){

    string prefixTrackFilename =  prefix_from_filename(filename);  
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");


    for (unsigned int i = 0; i < trackPts[k].size(); i++){
      //if (allTracks[k][i].valid == 1){
      if ((trackPts[k][i].valid == 1) && (trackPts[k][i].reflectance != 0) &&(trackPts[k][i].reflectance != -1)){//valid track and non-zero reflectance
        fprintf(fp, "%f\n", gain*trackPts[k][i].reflectance + bias);
      }
      else{
        fprintf(fp, "-1\n");
      }
    }

    fclose(fp);
    delete trackFilename;
  }

}

//saves the image points corresponding to a detectNum
void SaveImagePoints(vector< vector<LOLAShot> >  &allTracks, int detectNum, string filename)
{

  FILE *fp;

  for (unsigned int k = 0; k < allTracks.size(); k++ ){

    string prefixTrackFilename = prefix_from_filename(filename); 
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), k);
    fp = fopen(trackFilename, "w");

    for (unsigned int i = 0; i < allTracks[k].size()-1; i++){


      int found = 0;
      for (unsigned int j = 0; j < allTracks[k][i].LOLAPt.size(); j++){
        if ((allTracks[k][i].LOLAPt[j].s == detectNum) && (allTracks[k][i].valid == 1)){
          found = 1;
          fprintf(fp, "%f\n", allTracks[k][i].imgPt[j].val);
        }
      }
      if (found == 0){
        fprintf(fp, "-1\n");
      }

    }

    fclose(fp);
    delete trackFilename;
  }

}

//saves the altitude  points (LOLA) corresponding to a sensor ID = detectNum
void SaveAltitudePoints(vector< vector<LOLAShot> >  &tracks, int detectNum, string filename)
{

  FILE *fp;

  for (unsigned int t = 0; t < tracks.size(); t++ ){

    string prefixTrackFilename = prefix_from_filename(filename); 
    char* trackFilename = new char[500];
    sprintf (trackFilename, "%s_%d.txt", prefixTrackFilename.c_str(), t);
    fp = fopen(trackFilename, "w");

    for (unsigned int s = 0; s < tracks[t].size(); s++){

      int found = 0;
      for (unsigned int j = 0; j < tracks[t][s].LOLAPt.size(); j++){
         if ((tracks[t][s].LOLAPt[j].s == detectNum) && (tracks[t][s].valid == 1)){
             found = 1;
             fprintf(fp, "%f\n", tracks[t][s].LOLAPt[j].z());
         }
      }

      if (found == 0){
          fprintf(fp, "-1\n");
      }

    }

    fclose(fp);
    delete trackFilename;
  }

}


void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray, 
               string cubFile, vector<gcp> &gcpArray, Vector2 centroid)
{

   int index = 0;
    for (unsigned int t=0; t<trackPts.size(); t++){
      //int featureIndex = 0;
      for (unsigned int s=0; s<(unsigned int)trackPts[t].size(); s++){
	if ( (trackPts[t][s].featurePtLOLA == 1)/* && (trackPts[t][s].valid ==1 )*/){
          if (trackPts[t][s].valid == 1){          
	    gcpArray[index].filename.push_back(cubFile);
         
	    float i = (optimalTransfArray[0]*(trackPts[t][s].imgPt[2].x - centroid(0)) + 
                       optimalTransfArray[1]*(trackPts[t][s].imgPt[2].y - centroid(1)) + 
                       optimalTransfArray[2] + centroid(0));
	    float j = (optimalTransfArray[3]*(trackPts[t][s].imgPt[2].x - centroid(0)) + 
		       optimalTransfArray[4]*(trackPts[t][s].imgPt[2].y - centroid(1)) + 
                       optimalTransfArray[5] + centroid(1));
          
            gcpArray[index].x.push_back(i);
	    gcpArray[index].y.push_back(j);
	    gcpArray[index].x_before.push_back(trackPts[t][s].imgPt[2].x);
	    gcpArray[index].y_before.push_back(trackPts[t][s].imgPt[2].y);
            gcpArray[index].trackIndex = t;
            gcpArray[index].shotIndex = s;//featureIndex;
	  }
          //featureIndex++;
          index++;
	}
      }
    }
    
}

void UpdateGCP(vector<vector<LOLAShot> > trackPts, 
               vector<Vector4> matchArray, vector<float> errorArray, 
               string camCubFile, string mapCubFile, 
               vector<gcp> &gcpArray, float downsample_factor)
{

   std::vector<float> map_pixel;
   map_pixel.resize(2);
   
   std::vector<float> map_pixel_init;
   map_pixel_init.resize(2);

   std::vector<float> cam_pixel;
   cam_pixel.resize(2);
   
   std::vector<float> cam_pixel_init;
   cam_pixel_init.resize(2);

   Isis::Pvl label( mapCubFile);
   Isis::Camera* camera = Isis::CameraFactory::Create( label );
  // Note that ISIS is different from C. They start their index at 1.
  
   int featureIndex = 0;
   int validFeatureIndex = 0;

   for (unsigned int t = 0; t < trackPts.size(); t++){
     for (unsigned int s = 0; s < (unsigned int)trackPts[t].size(); s++){
       if (trackPts[t][s].featurePtLOLA == 1){
         if ((trackPts[t][s].valid == 1) && (trackPts[t][s].reflectance != 0) && (trackPts[t][s].reflectance != -1)){
          
	    gcpArray[featureIndex].filename.push_back(camCubFile);

            float i = matchArray[validFeatureIndex](2);
	    float j = matchArray[validFeatureIndex](3);
        
            map_pixel[0] = i;
            map_pixel[1] = j; 

            // convert a map projected pixel location to the
            // original image coordinate system.
             
            camera->SetImage(map_pixel[0], map_pixel[1]);
            cam_pixel[0] = camera->DetectorMap()->ParentSample();
            cam_pixel[1] = camera->DetectorMap()->ParentLine();
	   	
	    gcpArray[featureIndex].x.push_back(cam_pixel[0]/downsample_factor);
	    gcpArray[featureIndex].y.push_back(cam_pixel[1]/downsample_factor);

	    map_pixel_init[0]=trackPts[t][s].imgPt[2].x;
            map_pixel_init[1]=trackPts[t][s].imgPt[2].y;
            
            camera->SetImage(map_pixel_init[0], map_pixel_init[1]);
            cam_pixel_init[0] = camera->DetectorMap()->ParentSample();
            cam_pixel_init[1] = camera->DetectorMap()->ParentLine();
	       
            gcpArray[featureIndex].x_before.push_back(cam_pixel_init[0]/downsample_factor);
	    gcpArray[featureIndex].y_before.push_back(cam_pixel_init[1]/downsample_factor);
            gcpArray[featureIndex].trackIndex = t;
            gcpArray[featureIndex].shotIndex = s;

            validFeatureIndex++;
    
	    }//valid==1
	
	  featureIndex++;
	}
      }
    }

    // delete remaining ISIS objects
    delete camera;
    
}

void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray, 
               string camCubFile, string mapCubFile, vector<gcp> &gcpArray, Vector2 centroid, 
               float downsample_factor)
{

   std::vector<float> map_pixel;
   map_pixel.resize(2);
   
   std::vector<float> map_pixel_init;
   map_pixel_init.resize(2);

   std::vector<float> cam_pixel;
   cam_pixel.resize(2);
   
   std::vector<float> cam_pixel_init;
   cam_pixel_init.resize(2);

   Isis::Pvl label( mapCubFile);
   Isis::Camera* camera = Isis::CameraFactory::Create( label );
  // Note that ISIS is different from C. They start their index at 1.
  
   int index = 0;
   for (unsigned int t=0; t<trackPts.size(); t++){
     for (unsigned int s=0; s<(unsigned int)trackPts[t].size(); s++){
       if ((trackPts[t][s].featurePtLOLA == 1)/* && (trackPts[t][s].valid ==1)*/){
         if (trackPts[t][s].valid == 1){          
	    gcpArray[index].filename.push_back(camCubFile);
         
	    float i = (optimalTransfArray[0]*(trackPts[t][s].imgPt[2].x - centroid(0)) + 
                       optimalTransfArray[1]*(trackPts[t][s].imgPt[2].y - centroid(1)) + 
                       optimalTransfArray[2] + centroid(0));
	    float j = (optimalTransfArray[3]*(trackPts[t][s].imgPt[2].x - centroid(0)) + 
		       optimalTransfArray[4]*(trackPts[t][s].imgPt[2].y - centroid(1)) + 
                       optimalTransfArray[5] + centroid(1));

            map_pixel[0] = i;
            map_pixel[1] = j; 

            // convert a map projected pixel location to the
            // original image coordinate system.
             
            camera->SetImage(map_pixel[0], map_pixel[1]);
            cam_pixel[0] = camera->DetectorMap()->ParentSample();
            cam_pixel[1] = camera->DetectorMap()->ParentLine();
	   	
	    gcpArray[index].x.push_back(cam_pixel[0]/downsample_factor);
	    gcpArray[index].y.push_back(cam_pixel[1]/downsample_factor);

	    map_pixel_init[0]=trackPts[t][s].imgPt[2].x;
            map_pixel_init[1]=trackPts[t][s].imgPt[2].y;
            
            camera->SetImage(map_pixel_init[0], map_pixel_init[1]);
            cam_pixel_init[0] = camera->DetectorMap()->ParentSample();
            cam_pixel_init[1] = camera->DetectorMap()->ParentLine();
	       
            gcpArray[index].x_before.push_back(cam_pixel_init[0]/downsample_factor);
	    gcpArray[index].y_before.push_back(cam_pixel_init[1]/downsample_factor);
            gcpArray[index].trackIndex = t;
            gcpArray[index].shotIndex = s;
    
	  }//valid==1
	
	  index++;
	}
      }
    }

    // delete remaining ISIS objects
    delete camera;
    
}

void SaveGCPoints(vector<gcp> gcpArray,  string gcpFilename)
{

  for (unsigned int i = 0; i < gcpArray.size(); i++){
    
       //check if this GCP is valid
       int numFiles = gcpArray[i].filename.size();
       if (numFiles > 0){
	   stringstream ss;
	   ss<<i;

           stringstream featureIndexString;
           featureIndexString<<gcpArray[i].shotIndex;
           stringstream trackIndexString;
	   trackIndexString<<gcpArray[i].trackIndex;
           //string this_gcpFilename = gcpFilename+"_"+ss.str()+".gcp";
	   
           string this_gcpFilename = gcpFilename+"_"+trackIndexString.str()+"_"+featureIndexString.str()+".gcp";
	   FILE *fp = fopen(this_gcpFilename.c_str(), "w");
	 
	   fprintf(fp, "%f %f %f %f %f %f\n", 
		   gcpArray[i].lon, gcpArray[i].lat, gcpArray[i].rad, 
		   gcpArray[i].sigma_lon, gcpArray[i].sigma_lat, gcpArray[i].sigma_rad);
	 
	   for (int j = 0; j < numFiles-1; j++){
	     //string imgFilenameNoPath = sufix_from_filename(gcpArray[i].filename[j]);
              string imgFilenameNoPath = GetFilenameNoPath(gcpArray[i].filename[j]);
	      fprintf(fp,"%s %f %f\n", 
		      (char*)(imgFilenameNoPath.c_str()), gcpArray[i].x[j], gcpArray[i].y[j]);
	   }
	   if (numFiles > 0){
	     //string imgFilenameNoPath = sufix_from_filename(gcpArray[i].filename[numFiles-1]);
              string imgFilenameNoPath = GetFilenameNoPath(gcpArray[i].filename[numFiles-1]);
	      fprintf(fp, "%s %f %f", 
		      (char*)(imgFilenameNoPath.c_str()), gcpArray[i].x[numFiles-1], gcpArray[i].y[numFiles-1]);
	   }
	   fclose(fp);
       }
  }
   

}

//computes the min and max lon and lat of the LOLA data
//this function can be moved to a new util.c file
Vector4 FindMinMaxLat(vector<vector<LOLAShot> >trackPts)
{
  float minLat = 180;
  float maxLat = -180;
  float minLon = 180;
  float maxLon = -180;

  for (unsigned int i = 0; i < trackPts.size(); i++){
    for (unsigned int j = 0; j < trackPts[i].size(); j++){
      for(unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){
        float lon = trackPts[i][j].LOLAPt[k].x();
        float lat = trackPts[i][j].LOLAPt[k].y();

        if (lat < minLat){
          minLat = lat;
        }

        if (lon < minLon){
          minLon = lon;
        }

        if (lat > maxLat){
          maxLat = lat;
        }

        if (lon > maxLon){
          maxLon = lon;
        }
      }
    }
  }

  Vector4 coords;
  coords(0) = minLat; 
  coords(1) = maxLat; 
  coords(2) = minLon; 
  coords(3) = maxLon;

  return coords;

}

void ComputeAverageShotDistance(vector<vector<LOLAShot> >trackPts)
{
  int numValidPts = 0; 
  float avgDistance = 0.0;
  Datum moon;
  moon.set_well_known_datum("D_MOON");

  for (unsigned int i = 0; i < trackPts.size(); i++){
    for (unsigned int j = 1; j < trackPts[i].size(); j++){
      if ((trackPts[i][j].LOLAPt.size()==5) && (trackPts[i][j-1].LOLAPt.size()==5)){
         
        float lon, lat, rad;
        Vector3 lon_lat_rad;
        lon = trackPts[i][j-1].LOLAPt[0].x();
	lat = trackPts[i][j-1].LOLAPt[0].y();
	rad = trackPts[i][j-1].LOLAPt[0].z();
        lon_lat_rad(0) = lon;
        lon_lat_rad(1) = lat;
        //lon_lat_rad(2) = 1000*rad;
        //Vector3 prev_xyz = cartography::lon_lat_radius_to_xyz(lon_lat_rad);

        lon_lat_rad(2) = (rad-1737.4)*1000;
        Vector3 prev_xyz = moon.geodetic_to_cartesian(lon_lat_rad);
       
       
        lon = trackPts[i][j].LOLAPt[0].x();
	lat = trackPts[i][j].LOLAPt[0].y();
	rad = trackPts[i][j].LOLAPt[0].z();
        lon_lat_rad(0) = lon;
        lon_lat_rad(1) = lat;
        //lon_lat_rad(2) = 1000*rad;
        //Vector3 xyz = cartography::lon_lat_radius_to_xyz(lon_lat_rad);
        
        lon_lat_rad(2) = (rad-1737.4)*1000;
        Vector3 xyz = moon.geodetic_to_cartesian(lon_lat_rad);
        
	float dist_x = xyz[0]-prev_xyz[0];
        float dist_y = xyz[1]-prev_xyz[1];
        float dist_z = xyz[2]-prev_xyz[2];
        float dist = sqrt(dist_x*dist_x + dist_y*dist_y + dist_z*dist_z);
        avgDistance = avgDistance + dist;
        numValidPts = numValidPts + 1;

      }
    }
  }

  cout<<"numValidPts="<<numValidPts<<endl;  
  avgDistance = avgDistance/numValidPts;
  cout<<"avgDistance= "<<avgDistance<<endl;
}

//computes the average distance from the center of the shot to its neighbors in the shot
void ComputeAverageIntraShotDistance(vector<vector<LOLAShot> >trackPts)
{
  int numValidPts = 0; 
  // float avgDistance = 0.0;
  vector<Vector3> xyzArray;
  xyzArray.resize(5);
  vector<float> distArray;
  distArray.resize(5);

  for (unsigned int i = 0; i < trackPts.size(); i++){
    for (unsigned int j = 0; j < trackPts[i].size(); j++){
      if (trackPts[i][j].LOLAPt.size()==5){
        
        for (unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){
	  float lon, lat, rad;
	  Vector3 lon_lat_rad;
	  lon = trackPts[i][j].LOLAPt[k].x();
	  lat = trackPts[i][j].LOLAPt[k].y();
	  rad = trackPts[i][j].LOLAPt[k].z();
          //int s = trackPts[i][j].LOLAPt[0].s;

          //cout<<"s="<<s<<endl;

	  lon_lat_rad(0) = lon;
	  lon_lat_rad(1) = lat;
	  lon_lat_rad(2) = 1000*rad;
	  xyzArray[k] = cartography::lon_lat_radius_to_xyz(lon_lat_rad);
	}
     
        for (unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){
	  float x_dist = xyzArray[0][0]-xyzArray[k][0];
	  float y_dist = xyzArray[0][1]-xyzArray[k][1];
	  float z_dist = xyzArray[0][2]-xyzArray[k][2];
	  distArray[k] = distArray[k] + sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist);
	}

        numValidPts = numValidPts + 1;
        

      }
    }
  }

  cout<<"numValidPts="<<numValidPts<<endl;  
  for (unsigned int k = 0; k < distArray.size(); k++){
    distArray[k] = distArray[k]/numValidPts;
    cout<<k<<":distArray= "<<distArray[k]<<endl;
  }
 
}



int GetTimeDiff(pointCloud prevPt, pointCloud currPt, float timeThresh)
{


  double prevTime = prevPt.hour*3600.0 + prevPt.min*60.0 + prevPt.sec;
  double currTime = currPt.hour*3600.0 + currPt.min*60.0 + currPt.sec;
 
  if (prevPt.year != currPt.year){//new orbit
      return (1);
  }  

  if (prevPt.month != currPt.month){ //new orbit
      return (1);
  } 
  
  if (prevPt.day != currPt.day){//new orbit
      return (1);
  } 
 
  if (currTime - prevTime > timeThresh){ //new orbit
     return (1);
  }
  else{
     return (0);
  }
}
 



