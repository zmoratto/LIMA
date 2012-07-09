// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iomanip>
#include <math.h>
#include <boost/tokenizer.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <vw/Core.h>
#include <vw/Math/Matrix.h>
#include <vw/Cartography.h>
#include <vw/Cartography/SimplePointImageManipulation.h>
#include <vw/FileIO.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include <boost/filesystem.hpp>

#include "coregister.h"
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
pointCloud::pointCloud ( const Vector3& coords,
                         const int      y,
                         const int      mo,
                         const int      d,
                         const int      h,
                         const int      mi,
                         const float    se,
                         const int      detector 
                         ): Vector<double,3>( coords ) {
  year = y;
  month = mo;
  day = d;
  hour = h;
  min = mi;
  sec = se;
  s = detector;
}

void
LOLAShot::init
  (
  vector<pointCloud> pc,
  vector<imgPoint>   ip,
  vector<DEMPoint>   dp,
  int                v,
  // int                cpi,
  float              r,
  // float              si,
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
  // centerPtIndex = cpi;
  reflectance = r;
  // synthImage = si;
  calc_acp = ca;
  filter_response = fr;
  featurePtRefl = fpr;
  weightRefl = wr;
  featurePtLOLA = fpL;
  filresLOLA = fL;
  }

LOLAShot::LOLAShot( const vector<pointCloud>& pcv )
  {
  LOLAShot::init( pcv );
  }
LOLAShot::LOLAShot( const pointCloud& pt )
  {
  vector<pointCloud> pcv( 1, pt );
  LOLAShot::init( pcv );
  }

bool isTimeDiff( const pointCloud& p, const pointCloud& c, const float timeThresh){
  using namespace boost::posix_time;
  using namespace boost::gregorian;
  float intpart;
  ptime first( date(p.year,p.month,p.day),
               hours(p.hour) +minutes(p.min)+seconds(p.sec)
               +microseconds( modf(p.sec, &intpart) * 1000000 ) );
  ptime second( date(c.year,c.month,c.day),
                hours(c.hour)+minutes(c.min)+seconds(c.sec)
                +microseconds( modf(c.sec, &intpart) * 1000000 ) );

  if( first == second ) { return false; }

  time_duration duration( second - first );
  if( duration.is_negative() ) { duration = duration.invert_sign(); }

  if( duration.total_seconds() < timeThresh ){ return false; }
  return true;
}


vector<vector<LOLAShot> > LOLAFileRead( const string& f ) {
  // This is specifically for reading the RDR_*PointPerRow_csv_table.csv files only.
  // It will ignore lines which do not start with year (or a value that can be converted
  // into an integer greater than zero, specfically).
  char temp[255];
  ifstream file( f.c_str() );
  if( !file ) {
    vw_throw( vw::IOErr() << "Unable to open track file \"" << f << "\"" );
  }

  vector<vector<LOLAShot> > trackPts(1);

  // Prepare the tokenizer stuff.
  //typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  //boost::char_separator<char> sep(" ,");

  // file >> ignoreLine; // Skip header line
 
  // Establish 3d bounding box for acceptable LOLA data. 
  //Vector3 LOLAmin( -180, -90, 1720 );
  //Vector3 LOLAmax(  360,  90, 1750 );
  BBox3 LOLArange( Vector3(-180, -90, 1720), Vector3(360,  90, 1750) );

  int trackIndex = 0;
  pointCloud currPt;
  pointCloud prevPt;
  string line;
  while( getline(file, line, '\n') ) {

    // we went with C-style file reading instead of C++ in this instance
    // because we found it to be significantly faster on large files
    /*tokenizer tokens( line, sep );
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
    ++token;
    istringstream lon_s( *token );
    lon_s >> currPt.x(); // Pt_Longitude
    ++token;
    istringstream lat_s( *token );
    lat_s >> currPt.y(); // Pt_Latitude
    ++token;
    istringstream rad_s( *token );
    rad_s >> currPt.z(); // Pt_Radius
    advance( token, 8 );
    istringstream s_s( *token );
    s_s >>currPt.s; // S*/

    strncpy(temp, line.c_str(), 255);
    const char* token = strtok(temp, ",");
    
    int ret = sscanf(token, "%d-%d-%dT%d:%d:%g", &currPt.year, &currPt.month, &currPt.day, &currPt.hour, &currPt.min, &currPt.sec);
    if( currPt.year <= 0 ) { continue; }

    token = strtok(NULL, ",");
    ret += sscanf(token, "%lg", &currPt.x());

    token = strtok(NULL, ",");
    ret = ret + sscanf(token, "%lg", &currPt.y());
    token = strtok(NULL, ",");
    ret = ret + sscanf(token, "%lg", &currPt.z());

    for (int i = 0; i < 8; i++)
      token = strtok(NULL, ",");
    ret = ret + sscanf(token, "%d", &currPt.s);
    
    if (ret != 10)
    {
      fprintf(stderr, "Failed to read line %s.\n", line.c_str());
      vw_throw( vw::IOErr() << "Failed to read line." );
    }
    
    if( LOLArange.contains(currPt) ){
      if( trackPts[trackIndex].empty() ) {
        trackPts[trackIndex].push_back( LOLAShot(currPt) );
      }
      else {
        if( isTimeDiff(prevPt, currPt, 3000) ) { //new track
          trackIndex++;
          trackPts.resize(trackIndex+1); //start new track
          trackPts[trackIndex].push_back( LOLAShot(currPt) ); //put the new shot in it
        }
        else { //same track
          if( isTimeDiff(prevPt, currPt, 0) ) { //new shot
            trackPts[trackIndex].push_back( LOLAShot(currPt) );
          }
          else { //same shot
            trackPts[trackIndex].back().LOLAPt.push_back( currPt );
          }
        }
      }

      //copy current pc into prevPt
      prevPt = currPt;
    }
  } 
  file.close();
  return trackPts; 
}

vector<vector<LOLAShot> > LOLAFilesRead( const string& tracksList )
{
	vector<string> trackFiles = ReadVectorFrom<string>(tracksList);
	vector<vector<LOLAShot> > points;
	for (unsigned int i = 0; i < trackFiles.size(); i++)
	{
		vector<vector<LOLAShot> > shots = LOLAFileRead(trackFiles[i]);
		points.reserve(points.size() + distance(shots.begin(),shots.end()));
		points.insert(points.end(), shots.begin(), shots.end());
	}

	return points;
}

Vector2 ComputeMinMaxValuesFromCub( ImageView<PixelGray<float> > & isis_view)
{
  Vector2 minmax;
  PixelGray<float> minVal = std::numeric_limits<float>::max();
  PixelGray<float> maxVal = std::numeric_limits<float>::min();
  min_max_pixel_values( isis_view, minVal, maxVal );
  minmax(0) = minVal;
  minmax(1) = maxVal;
  //cout<<"min="<<minVal<<", max="<<maxVal<<endl;
  return minmax;
}

/*
Vector2 ComputeMinMaxValuesFromDEM(string demFilename)
{
  Vector2 minmax;
  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(demFilename) );
  double nodataVal = rsrc->nodata_read();
  //cout<<"nodaval:"<<nodataVal<<endl;
  DiskImageView<PixelGray<uint16> > dem( rsrc );
  int width = dem.cols();
  int height = dem.rows();
  float minVal = 100000000.0;
  float maxVal = -100000000.0;
  for (int i = 0; i < height; i++){
    for (int j = 0; j < width; j++){
      
      if ((dem(j,i) < minVal) && (dem(j,i) > nodataVal)){
	minVal = dem(j,i);
      }
      if ((dem(j,i) > maxVal) && (dem(j,i) > nodataVal)){
	maxVal = dem(j,i);
      }
    }
  }
  minmax(0) = minVal;
  minmax(1) = maxVal;
  cout<<"min="<<minVal<<", max="<<maxVal<<endl;

  return minmax;
}
*/
int GetAllPtsFromCub( vector<vector<LOLAShot> >& trackPts, camera::IsisCameraModel & model, ImageView<PixelGray<float> > & cubImage)
{
  //calculate the min max value of the image
  Vector2 minmax = ComputeMinMaxValuesFromCub( cubImage );

  // float minTrackVal = std::numeric_limits<float>::max();
  // float maxTrackVal = std::numeric_limits<float>::min();

  int numValidImgPts = 0;

  InterpolationView<EdgeExtensionView<ImageView<PixelGray<float> >, ConstantEdgeExtension>, BilinearInterpolation> interpImg
    = interpolate(cubImage, BilinearInterpolation(), ConstantEdgeExtension());

  BBox2i bbox = bounding_box( cubImage );  
  for( unsigned int k = 0; k < trackPts.size(); ++k ){
    for( unsigned int i = 0; i < trackPts[k].size(); ++i ){
      vector<pointCloud> points = trackPts[k][i].LOLAPt;
      
      trackPts[k][i].valid = 1; 
      trackPts[k][i].imgPt.resize( points.size() );

      for( unsigned int j = 0; j < points.size(); ++j ){
        Vector3 lon_lat_rad( points[j].x(), points[j].y(), points[j].z()*1000 );
        Vector3 xyz = cartography::lon_lat_radius_to_xyz(lon_lat_rad);
        Vector2 cub_pix = model.point_to_pixel(xyz);
        //check that (x,y) are within the image boundaries and for validity
        if( (bbox.contains(cub_pix)) && 
            (interpImg( cub_pix.x(), cub_pix.y() ) > minmax(0)) ){//valid values
          trackPts[k][i].imgPt[j].val = interpImg( cub_pix.x(), cub_pix.y() );
          trackPts[k][i].imgPt[j].x = cub_pix.x();
          trackPts[k][i].imgPt[j].y = cub_pix.y();
          // if( interpImg(x,y) < minTrackVal ){ minTrackVal = interpImg(x,y); }
          // if( interpImg(x,y) > maxTrackVal ){ maxTrackVal = interpImg(x,y); }
          numValidImgPts++;
        }
        else{//invalidate the point
          trackPts[k][i].valid = 0;
          break;
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
  //cout <<"minTrackVal="<<minTrackVal<<", maxTrackVal="<<maxTrackVal<<endl;
  
  return numValidImgPts;
}


/* Deprecated?  This is commented out in the .h file.
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
            nominator = nominator + (trackPts[k][i].imgPt[j].val//-minmax(0)
                                    )/trackPts[k][i].reflectance;
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
*/

void GainBiasAccumulator( const vector<AlignedLOLAShot>& trackPts, 
                                float&            sum_rfl, 
                                float&            sum_img, 
                                float&            sum_rfl_2, 
                                float&            sum_rfl_img, 
                                int&              numValidPts ){
  for (unsigned int i = 0; i < trackPts.size(); ++i){
    if ((trackPts[i].reflectance >= 0) && trackPts[i].image >= 0){
      //update the nominator for the center point
	      sum_rfl     += trackPts[i].reflectance;
	      sum_rfl_2   += trackPts[i].reflectance*trackPts[i].reflectance;
	      sum_rfl_img += trackPts[i].reflectance*(trackPts[i].image);
	      sum_img     += trackPts[i].image;
	      ++numValidPts;
      }
    }
}

Vector2 GainBiasSolver( const float& sum_rfl, 
                        const float& sum_img, 
                        const float& sum_rfl_2, 
                        const float& sum_rfl_img, 
                        const int&   numValidPts ){
  /*Matrix<float,2,2> rhs;
  Vector<float,2> lhs;

  rhs(0,0) = sum_rfl_2;
  rhs(0,1) = sum_rfl;
  rhs(1,0) = sum_rfl;
  rhs(1,1) = numValidPts;
  lhs(0) = sum_rfl_img;
  lhs(1) = sum_img;
  solve_symmetric_nocopy(rhs,lhs);
  return lhs;*/
  Vector<float,2> result;
  result(0) = sum_img / sum_rfl;
  result(1) = 0;
  return result;
}


//computes the gain and bias factor for each track
Vector2 ComputeGainBiasFactor( const vector<AlignedLOLAShot>& trackPts ) {
  int numValidPts = 0;
  float sum_rfl = 0.0; 
  float sum_img = 0.0;
  float sum_rfl_2 = 0.0;
  float sum_rfl_img = 0.0;
  Vector2 gain_bias;

  GainBiasAccumulator( trackPts, sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );

  //cout<<"NUM_VALID_POINTS="<<numValidPts<<endl;
  //if (numValidPts != 0){
  if( numValidPts > 1 ){ 
    gain_bias = GainBiasSolver( sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
  }
  else{
    //invalid scaleFactor, all tracks are invalid
    gain_bias(0) = 0.0;
    gain_bias(1) = 0.0;
  }
  return gain_bias;
}

//computes the gain and bias factor for all tracks at once
Vector2 ComputeGainBiasFactor( const vector<vector<AlignedLOLAShot> >& trackPts ) {
  int numValidPts = 0;
  float sum_rfl = 0.0; 
  float sum_img = 0.0;
  float sum_rfl_2 = 0.0;
  float sum_rfl_img = 0.0;
  Vector2 gain_bias;

  for( unsigned int k = 0; k < trackPts.size(); ++k){
    GainBiasAccumulator( trackPts[k], sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
  }

  //cout<<"NUM_VALID_POINTS="<<numValidPts<<endl;
  if (numValidPts != 0){ 
    gain_bias = GainBiasSolver( sum_rfl, sum_img, sum_rfl_2, sum_rfl_img, numValidPts );
  }
  else{
    //invalid scaleFactor, all tracks are invalid
    gain_bias(0) = 0.0;
    gain_bias(1) = 0.0;
  }
  return gain_bias;
}

/* Deprecated?  Not used anywhere
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
*/

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

int ComputeAllReflectance(       vector< vector<LOLAShot> >& shots,  
                           const Vector3&                    cameraPosition, 
                           const Vector3&                    lightPosition) {
  int numValidReflPts = 0;

  float minReflectance = std::numeric_limits<float>::max();
  float maxReflectance = std::numeric_limits<float>::min();
 
  for( unsigned int k = 0; k < shots.size(); ++k ){
    for( unsigned int i = 0; i < shots[k].size(); ++i ){
      try {
        if( shots[k][i].LOLAPt.size() > 5 ) { 
          vw_throw( ArgumentErr() << "Too many LOLAPts: " << shots[k][i].LOLAPt.size() );
        }
        pointCloud centerPt = GetPointFromIndex( shots[k][i].LOLAPt, 1 );
        pointCloud topPt    = GetPointFromIndex( shots[k][i].LOLAPt, 3 );
        pointCloud leftPt   = GetPointFromIndex( shots[k][i].LOLAPt, 2 );
      
        centerPt.z() *= 1000;
        topPt.z()    *= 1000;
        leftPt.z()   *= 1000;
        
        Vector3 xyz = lon_lat_radius_to_xyz(centerPt);
        Vector3 xyzTop = lon_lat_radius_to_xyz(topPt);
        Vector3 xyzLeft = lon_lat_radius_to_xyz(leftPt);

        Vector3 normal = ComputeNormalFrom3DPointsGeneral(xyz, xyzLeft, xyzTop);

        shots[k][i].reflectance = ComputeLunarLambertianReflectanceFromNormal(lightPosition, cameraPosition, xyz, normal); 
        
        if( shots[k][i].reflectance != 0 ){ 
          ++numValidReflPts;

          if( shots[k][i].reflectance < minReflectance ){
            minReflectance = shots[k][i].reflectance;
          } 
          if( shots[k][i].reflectance > maxReflectance ){
            maxReflectance = shots[k][i].reflectance;
          }
        }        
      }
      catch( const vw::ArgumentErr& error ){ shots[k][i].reflectance = -1; }
    }//i
  }//k
 
  //cout<<"NEW:"<<"minReflectance="<<minReflectance<<", maxReflectance="<<maxReflectance<<endl;

  return numValidReflPts;
}

pointCloud GetPointFromIndex( const vector<pointCloud>& LOLAPts, const int index ) {
  for( unsigned int i = 0; i < LOLAPts.size(); ++i ) {
    if( LOLAPts[i].s == index ){
      return LOLAPts[i];
    }
  }
  vw_throw( ArgumentErr() << "Couldn't find a point with detector number " << index );
}


void SaveReflectancePoints( const vector<vector<LOLAShot> >& trackPts, 
                            const Vector2&                   gain_bias,
                            const string&                    filename) {
  boost::filesystem::path p( filename );
  for( unsigned int k = 0; k < trackPts.size(); ++k ){
    ostringstream os;
    os << p.stem().string() << "_" << k << ".txt";

    ofstream file( os.str().c_str() );
    if( !file ) {
      vw_throw( ArgumentErr() << "Can't open reflectance output file " << os.str() );
    }
    
    for( unsigned int i = 0; i < trackPts[k].size(); ++i){
      if ( (trackPts[k][i].valid == 1) && 
           (trackPts[k][i].reflectance > 0) ){//valid track and non-zero reflectance
        file << ( gain_bias(0)*trackPts[k][i].reflectance + gain_bias(1) ) << endl;
      }
      else{ file << "-1" << endl; }
    }
    file.close();
  }
}

//saves the image points corresponding to a detectNum
void SaveImagePoints( const vector<vector<LOLAShot> >& allTracks,
                      const int                        detectNum,
                      const string&                    filename) {
  boost::filesystem::path p( filename );
  for (unsigned int k = 0; k < allTracks.size(); ++k ){
    ostringstream os;
    os << p.stem().string() << "_" << k << ".txt";
    
    ofstream file( os.str().c_str() );
    if( !file ) {
      vw_throw( ArgumentErr() << "Can't open image point output file " << os.str() );
    }

    for( unsigned int i = 0; i < allTracks[k].size(); ++i){
      bool found = false;
      for( unsigned int j = 0; j < allTracks[k][i].LOLAPt.size(); ++j){
        if( (allTracks[k][i].LOLAPt[j].s == detectNum) && 
            (allTracks[k][i].valid       == 1)           ){
          found = true;
          file <<  allTracks[k][i].imgPt[j].val << endl;
        }
      }
      if( !found ){ file << "-1" << endl; }
    }
    file.close();
  }
}

//saves the altitude  points (LOLA) corresponding to a sensor ID = detectNum
void SaveAltitudePoints( const vector<vector<LOLAShot> >& tracks,
                         const int                        detectNum,
                         const string&                    filename) {
  boost::filesystem::path p( filename );
  for (unsigned int t = 0; t < tracks.size(); ++t ){
    ostringstream os;
    os << p.stem().string() << "_" << t << ".txt";
    
    ofstream file( os.str().c_str() );
    if( !file ) {
      vw_throw( ArgumentErr() << "Can't open altitude output file " << os.str() );
    }

    for( unsigned int s = 0; s < tracks[t].size(); ++s ){
      bool found = false;
      for( unsigned int j = 0; j < tracks[t][s].LOLAPt.size(); ++j ){
        if( (tracks[t][s].LOLAPt[j].s == detectNum) && 
            (tracks[t][s].valid       == 1)           ){
          found = true;
          file << tracks[t][s].LOLAPt[j].z() << endl;
        }
      }
      if( !found ){ file << "-1" << endl; }
    }
    file.close();
  }
}

/* Deprecated?
void UpdateGCP(vector<vector<LOLAShot> > trackPts, Vector<float, 6> optimalTransfArray, 
               string cubFile, vector<gcp> &gcpArray, Vector2 centroid)
{

   int index = 0;
    for (unsigned int t=0; t<trackPts.size(); t++){
      //int featureIndex = 0;
      for (unsigned int s=0; s<(unsigned int)trackPts[t].size(); s++){
	if ( (trackPts[t][s].featurePtLOLA == 1)
          // && (trackPts[t][s].valid ==1 )
          ){
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
*/

// update ground control points
void UpdateGCP( const vector<vector<LOLAShot> >& trackPts, 
                const vector<Vector4>&           matchArray, 
                const string&                    camCubFile, 
                const string&                    mapCubFile, 
                vector<gcp>&                     gcpArray,
		float downsample_factor) {
  int featureIndex = 0;
  int validFeatureIndex = 0;
  Isis::Pvl label( mapCubFile);
  Isis::Camera* camera = Isis::CameraFactory::Create( label );

  for( unsigned int t = 0; t < trackPts.size(); ++t ){
    for( unsigned int s = 0; s < trackPts[t].size(); ++s ){
      if( trackPts[t][s].featurePtLOLA == 1 ){
        gcpArray[featureIndex].filename.push_back(camCubFile);
        if( (trackPts[t][s].valid == 1) && 
            (trackPts[t][s].reflectance != 0) && 
            (trackPts[t][s].reflectance != -1)  ){

          // convert a map projected pixel location to the
          // original image coordinate system.
          camera->SetImage( matchArray[validFeatureIndex](2), matchArray[validFeatureIndex](3) );
          gcpArray[featureIndex].x.push_back( ( camera->DetectorMap()->ParentSample() ) / downsample_factor );
          gcpArray[featureIndex].y.push_back( ( camera->DetectorMap()->ParentLine() ) / downsample_factor );

          camera->SetImage( trackPts[t][s].imgPt[2].x, trackPts[t][s].imgPt[2].y );
          gcpArray[featureIndex].x_before.push_back(
                                 ( camera->DetectorMap()->ParentSample() )/downsample_factor );
          gcpArray[featureIndex].y_before.push_back(
                                 ( camera->DetectorMap()->ParentLine() )/downsample_factor );

          gcpArray[featureIndex].trackIndex = t;
          gcpArray[featureIndex].shotIndex = s;

          validFeatureIndex++;
        }//valid==1
	else
	{
          	gcpArray[featureIndex].x.push_back(-1);
          	gcpArray[featureIndex].y.push_back(-1);
          	gcpArray[featureIndex].x_before.push_back(-1);
          	gcpArray[featureIndex].y_before.push_back(-1);
		//gcpArray[featureIndex].trackIndex = -1;
		//gcpArray[featureIndex].shotIndex = -1;
	}
        featureIndex++;
      }
    }
  }
}

/* Deprecated?
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
       if ((trackPts[t][s].featurePtLOLA == 1)
         // && (trackPts[t][s].valid ==1)
         ){
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
*/

void SaveGCPoints( const vector<gcp>& gcpArray, const string& gcpFilename ) {
  for( unsigned int i = 0; i < gcpArray.size(); ++i ){
    //check if this GCP is valid
    if( !gcpArray[i].filename.empty() ){
      stringstream filename;
      filename << gcpFilename << "_" << gcpArray[i].trackIndex 
                              << "_" << gcpArray[i].shotIndex << ".gcp";
      
      ofstream file( filename.str().c_str() );
      if( !file ) {
        vw_throw( ArgumentErr() << "Can't open gcp output file " << filename.str() );
      }

      file << fixed << setprecision(6)
           << gcpArray[i].lon << " "
           << gcpArray[i].lat << " "
           << gcpArray[i].rad << " "
           << gcpArray[i].sigma_lon << " "
           << gcpArray[i].sigma_lat << " "
           << gcpArray[i].sigma_rad << " " << endl;

      for( unsigned int j = 0; j < gcpArray[i].filename.size(); ++j ){
        boost::filesystem::path p( gcpArray[i].filename[j] );
        file << p.filename().string() << " " 
             << fixed << setprecision(6) << gcpArray[i].x[j] << " " << gcpArray[i].y[j] << endl;
	  }
	  file.close();
    }
  }
}

//computes the min and max lon and lat of the LOLA data
BBox2 FindShotBounds(const vector<vector<LOLAShot> >& trackPts) {
  Vector<float,2> min( std::numeric_limits<float>::max(), 
                       std::numeric_limits<float>::max() );
  Vector<float,2> max( std::numeric_limits<float>::min(), 
                       std::numeric_limits<float>::min() );

  for(     unsigned int i = 0; i < trackPts.size();              ++i){
    for(   unsigned int j = 0; j < trackPts[i].size();           ++j){
      for( unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); ++k){
        float lon = trackPts[i][j].LOLAPt[k].x();
        float lat = trackPts[i][j].LOLAPt[k].y();

        if (lat < min.y() ){ min.y() = lat; }
        if (lon < min.x() ){ min.x() = lon; }
        if (lat > max.y() ){ max.y() = lat; }
        if (lon > max.x() ){ max.x() = lon; }
      }
    }
  }

  BBox2 bounds( min, max );
  return bounds;
}

Vector4 FindMinMaxLat( const vector<vector<LOLAShot> >& trackPts ) {
  BBox2 bounds = FindShotBounds( trackPts );
  Vector4 coords;
  coords(0) = bounds.min().y();
  coords(1) = bounds.max().y();
  coords(2) = bounds.min().x();
  coords(3) = bounds.max().x();
  return coords;
}

/* This is not currently being used, but may be useful in the future.
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
*/

/* This is not currently being used, but may be useful in the future.
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
*/

void apply_gain_bias(vector<vector<AlignedLOLAShot> >& tracks)
{
	Vector2 gain_bias = ComputeGainBiasFactor( tracks );
	for (unsigned int i = 0; i < tracks.size(); i++)
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			tracks[i][j].synth_image = tracks[i][j].reflectance * gain_bias(0) + gain_bias(1);
			if (tracks[i][j].synth_image > 1.0)
				tracks[i][j].synth_image = 1.0;
		}
}

void apply_gain_bias(vector<AlignedLOLAShot>& track)
{
	Vector2 gain_bias = ComputeGainBiasFactor( track );
	for (unsigned int i = 0; i < track.size(); i++)
	{
		track[i].synth_image = track[i].reflectance * gain_bias(0) + gain_bias(1);
		if (track[i].synth_image > 1.0)
			track[i].synth_image = 1.0;
	}
}

void transform_track_no_gain_bias(vector<AlignedLOLAShot> & track, Matrix3x3 transform, ImageView<PixelGray<float> >& cub)
{
	for (unsigned int i = 0; i < track.size(); i++)
	{
		if (track[i].reflectance == -1 || track[i].imgPt.size() <= 2)
		{
			track[i].image = -1;
			continue;
		}
		Vector3 r = transform * Vector3(track[i].imgPt[2].x, track[i].imgPt[2].y, 1);
		int x = (int)r(0), y = (int)r(1);
		track[i].image_x = x;
		track[i].image_y = y;
		if (x < 0 || y < 0 || x >= cub.cols() || y >= cub.rows())
		{
			track[i].image = -1;
			continue;
		}
		track[i].image = cub(x, y)[0];
	}
}

void transform_track(vector<AlignedLOLAShot> & track, Matrix3x3 transform, ImageView<PixelGray<float> >& cub)
{
	transform_track_no_gain_bias(track, transform, cub);
	
	apply_gain_bias(track);
}


void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Matrix3x3 transform, ImageView<PixelGray<float> >& cub)
{
	for (unsigned int i = 0; i < tracks.size(); i++)
		transform_track_no_gain_bias(tracks[i], transform, cub);

	apply_gain_bias(tracks);
}

void transform_tracks(vector<vector<AlignedLOLAShot> > & tracks, Matrix3x3 transform, string cubFile)
{
	boost::shared_ptr<DiskImageResource> rsrc(new DiskImageResourceIsis(cubFile));
	DiskImageView<PixelGray<float> > cub(rsrc);
	ImageView<PixelGray<float> > assembledImg(cub.cols(), cub.rows());
	assembledImg = normalize(cub);
	transform_tracks(tracks, transform, assembledImg);
}

// this shifts the *original* image points
void transform_tracks_by_matrices(vector<vector<LOLAShot> > & tracks, vector<Matrix3x3> matrices)
{
	for (unsigned int i = 0; i < tracks.size(); i++)
	{
		Matrix3x3 M = matrices[i];
		for (unsigned int j = 0; j < tracks[i].size(); j++)
			for (unsigned int k = 0; k < tracks[i][j].imgPt.size(); k++)
			{
				Vector3 r = M * Vector3(tracks[i][j].imgPt[k].x, tracks[i][j].imgPt[k].y, 1);
				tracks[i][j].imgPt[k].x = r(0);
				tracks[i][j].imgPt[k].y = r(1);
			}
	}
}

void save_track_data( const std::vector<std::vector<AlignedLOLAShot> >& tracks, const std::string& filename)
{
	FILE* f = fopen(filename.c_str(), "w");
	if (f == NULL)
	{
		fprintf(stderr, "Could not open track data file %s for writing.\n", filename.c_str());
		return;
	}
	for (unsigned int i = 0; i < tracks.size(); i++)
	{
		int num_valid = 0;
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0)
				continue;
			num_valid++;
		}
		// ignore small tracks
		if (num_valid < 25)
			continue;
		// first print altitude, then distance, then synthetic image, then corresponding image pixel
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0)
				continue;
			fprintf(f, "%g", tracks[i][j].LOLAPt[2][2]);
			if (j != tracks[i].size() - 1)
				fprintf(f, " ");
		}
		fprintf(f, "\n");
		float last_x = 0.0, last_y = 0.0;
		float total_distance = 0.0;
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0)
				continue;
        		Vector3 lon_lat_rad( tracks[i][j].LOLAPt[2].x(), tracks[i][j].LOLAPt[2].y(), tracks[i][j].LOLAPt[2].z()*1000 );
        		Vector3 xyz = cartography::lon_lat_radius_to_xyz(lon_lat_rad);
			if (j != 0)
				total_distance += sqrt(pow(xyz.x() - last_x, 2) + pow(xyz.y() - last_y, 2));
			fprintf(f, "%g", total_distance);
			if (j != tracks[i].size() - 1)
				fprintf(f, " ");
			last_x = xyz.x();
			last_y = xyz.y();
		}
		fprintf(f, "\n");
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0)
				continue;
			fprintf(f, "%g", tracks[i][j].synth_image);
			if (j != tracks[i].size() - 1)
				fprintf(f, " ");
		}
		fprintf(f, "\n");
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			if (tracks[i][j].synth_image < 0 || tracks[i][j].image < 0)
				continue;
			fprintf(f, "%g", tracks[i][j].image);
			if (j != tracks[i].size() - 1)
				fprintf(f, " ");
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

vector<Matrix3x3> load_track_transforms(const std::string& filename)
{
	vector<Matrix3x3> l;
	FILE* f = fopen(filename.c_str(), "r");

	while (true)
	{
		Matrix3x3 m((double)1.0, 0, 0, 0, 1, 0, 0, 0, 1);
		float t[6];
		int ret = fscanf(f, "%g %g %g %g %g %g\n", &t[0], &t[1], &t[2], &t[3], &t[4], &t[5]);
		if (ret != 6)
			return l;
		m(0, 0) = t[0]; m(0, 1) = t[1]; m(0, 2) = t[2];
		m(1, 0) = t[3]; m(1, 1) = t[4]; m(1, 2) = t[5];
		l.push_back(m);
	}

	return l;
}

std::vector<std::vector< AlignedLOLAShot> > initialize_aligned_lola_shots(std::vector<std::vector<LOLAShot> >& trackPts)
{
	std::vector<std::vector< AlignedLOLAShot> > tracks;
	for (unsigned int i = 0; i < trackPts.size(); i++)
	{
		std::vector< AlignedLOLAShot> t;
		for (unsigned int j = 0; j < trackPts[i].size(); j++)
		{
			if (!trackPts[i][j].valid)
				continue;
			AlignedLOLAShot s(trackPts[i][j]);
			if (trackPts[i][j].imgPt.size() > 2)
			{
				s.image_x = trackPts[i][j].imgPt[2].x;
				s.image_y = trackPts[i][j].imgPt[2].y;
				s.image = trackPts[i][j].imgPt[2].val;
			}
			t.push_back(s);
		}
		tracks.push_back(t);
	}
	return tracks;
}
