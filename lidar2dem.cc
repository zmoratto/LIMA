// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

//boost
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

//VW
#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Math/Matrix.h>

//ATK
#include "tracks.h"
#include "coregister.h"
#include "icp.h"
#include "util.h"

using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

int main( int argc, char *argv[] )
{
 int verbose;
 std::string inputCSVFilename; 
 std::string configFilename;
 std::string errorFilename;
 std::vector<std::string> inputDEMFiles;
 std::string precDEMDir="NO_DIR";
 std::string resDir;
 std::string auxDir;

 po::options_description general_options("Options");
 general_options.add_options()
    ("Lidar-filename,l", po::value<std::string>(&inputCSVFilename))
    ("inputDEMFiles,d", po::value<std::vector<std::string> >(&inputDEMFiles))
    ("DEM-precision-directory,p", 
     po::value<std::string>(&precDEMDir), 
		"DEM precision directory.") 
    ("results-directory,r", 
		po::value<std::string>(&resDir)->default_value("../results"), 
		"results directory.") // Currently a no-op until we implement it below.
    ("settings-filename,s", 
		po::value<std::string>(&configFilename)->default_value("lidar2dem_settings.txt"), 
		"settings filename.")
    ("error-filename,e", 
		po::value<std::string>(&errorFilename)->default_value("lidar2dem_errors.txt"), 
		"error output filename.")
    ("verbose,v", 
		po::value<int>(&verbose)->default_value(1), 
		"Verbosity level, zero emits no messages.")
    ("help,h", "Display this help message");
  

  po::options_description hidden_options("");

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("inputDEMFiles", -1);

  std::ostringstream usage;
  usage << "Description: main code for Lidar to DEM co-registration" << std::endl << std::endl;
  usage << general_options << std::endl;
 
  po::variables_map vm;
  try
    {
      po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
      po::notify( vm );
    }
  catch ( po::error &e )
   {
     std::cout << "An error occured while parsing command line arguments.\n";
     std::cout << "\t" << e.what() << "\n\n";
     std::cout << usage.str();
     return 1;
   }
  
  if( vm.count("help") )
    {
      std::cerr << usage.str() << std::endl;
      return 1;
    }

  if (vm.count(inputCSVFilename)){
      std::cerr << "Error: Must specify at least  one Lidar file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
  }

  vector<string> DEMFiles; 
  int numDEMInputFiles = inputDEMFiles.size();

  if(numDEMInputFiles < 1) {
    std::cerr << "Error: Must specify at least  one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }

  //determine if inputDEMFiles is a text file containing a list of DEM files, one DEM file or a set of DEM files
  //by reading the file extension. A text file containing the DEM list *must* have extension .txt 
  DEMFiles = AccessDataFilesFromInput(inputDEMFiles);
 
  //Set up VW logging
  vw_log().console_log().rule_set().add_rule(vw::InfoMessage,"*");

 
  struct CoregistrationParams settings;
  if( ReadConfigFile(configFilename, &settings) && verbose > 0 ){
   cout << "Config file " << configFilename << " found." << endl;
  } 
  else if( verbose > 0 ){
   cout << "Config file " << configFilename << " not found, using defaults." << endl;
  }

  if( verbose > 0 ){ 
   cout << settings << endl; 
  }


  int numDEMFiles = DEMFiles.size();
  cout<<"numDEMFiles="<<numDEMFiles<<endl;

  auxDir = resDir+"/aux";
  cout<<"aux_dir="<<auxDir<<endl;
  //create the results and auxiliary directories
  string makeResDirCmd = "mkdir " + resDir;
  string makeAuxDirCmd = "mkdir " + auxDir;
  system(makeResDirCmd.c_str()); 
  system(makeAuxDirCmd.c_str());


  //read LOLA tracks
  vector<vector<LOLAShot> > trackPts;
  try{
    trackPts =  CSVFileRead(inputCSVFilename);
    if( verbose > 0 ){ cout << "Tracks file: " << inputCSVFilename << endl; }
  }
  catch (const vw::IOErr& error){
    cerr << error.what() << endl;
    exit(1);
  }
 
  //determine the overlapping DEMs - START
  cout<<"Selecting the overlapping DEMs ..."<<endl;
  std::vector<int> overlapIndices;
  string overlapListFilename = auxDir+string("/")+GetFilenameNoPath(inputCSVFilename);
  FindAndReplace(overlapListFilename, ".csv", "_overlap_list.txt"); 
  int fileFound = ReadOverlapList(overlapListFilename, overlapIndices);
  if (fileFound == 0){
    Vector4 lat_lon_bb = FindMinMaxLat(trackPts); 
    Vector4 lon_lat_bb;
    lon_lat_bb[0]=lat_lon_bb[2];
    lon_lat_bb[1]=lat_lon_bb[3];
    lon_lat_bb[2]=lat_lon_bb[0];
    lon_lat_bb[3]=lat_lon_bb[1];
    //printf("lidar corners: %f %f %f %f\n", lon_lat_bb[0], lon_lat_bb[1], lon_lat_bb[2], lon_lat_bb[3]);
    overlapIndices = makeOverlapListFromGeoTiff(DEMFiles, lon_lat_bb);
    SaveOverlapList(overlapListFilename, overlapIndices);
  }
  PrintOverlapList(overlapIndices);
  int numOverlappingDEMs;
  if ((overlapIndices.size() == 1) && (overlapIndices[0] == -1)){
    numOverlappingDEMs = 0;
  }
  else{
     numOverlappingDEMs = (int)overlapIndices.size();
  }
  cout<<"numOverlapDEMs="<<numOverlappingDEMs<<endl;
  cout<<"done."<<endl;
  //determine the overlapping DEMs - END
   
  for (int index = 0; index < numOverlappingDEMs; index++){
    
    string inputDEMFilename = DEMFiles[overlapIndices[index]];

    if (verbose > 0 ){ 
      cout <<"DEM filename: " << inputDEMFilename << endl; 
    }
    
    //read DEM file      
    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL(inputDEMFilename) );
    if (rsrc->has_nodata_read()){
        settings.noDataVal = rsrc->nodata_read();
	if( verbose > 0 ){
	  cout 	<< "Found nodata value, " << settings.noDataVal 
		<< ", in " << inputDEMFilename << endl;
	}
    }
    else if( verbose > 0 ){
      cout << "Using default nodata value: " << settings.noDataVal << endl;
    }
    
    DiskImageView<float> DEM(rsrc);
    GeoReference DEMGeo;
    read_georeference(DEMGeo, inputDEMFilename);

    if (precDEMDir.compare(string("NO_DIR"))!=0){
      //select DEM points with high precision closest to LOLA tracks
      string inputPrecFilename = GetFilenameNoPath(inputDEMFilename);
      FindAndReplace(inputPrecFilename, string("dem"), string("Prec"));
      inputPrecFilename = precDEMDir + string("/") + inputPrecFilename;
      cout<<"DEMFilename="<<inputDEMFilename<<endl;
      cout<<"PrecFilename="<<inputPrecFilename<<endl;
      boost::shared_ptr<DiskImageResource> prec_rsrc( new DiskImageResourceGDAL(inputPrecFilename) );
      DiskImageView<float> DEM_Prec(prec_rsrc);
      //select DEM points closest to LOLA tracks
      GetAllPtsFromDEM_Prec(trackPts, DEM, DEMGeo, settings.noDataVal, DEM_Prec);
    }  
    else{
      //select DEM points closest to LOLA tracks
      GetAllPtsFromDEM(trackPts, DEM, DEMGeo, settings.noDataVal);
    }

    Vector3 currTranslation;
    Matrix<float, 3,3 > currRotation;
    vector<Vector3> xyzModelArray;//LOLA
    vector<Vector3> llrModelArray;//LOLA
    vector<Vector3> xyzMatchArray;//DEM
    valarray<float> xyzErrorArray;//error array same size as model and feature
    vector<float>   radErrorArray; //altitude-radial error array
    
    currRotation[0][0] = 1.0;
    currRotation[1][1] = 1.0;
    currRotation[2][2] = 1.0;
    
    //copy info to featureArray and modelArray
    for(unsigned int k = 0; k < trackPts.size();k++){
      for(unsigned int i = 0; i < trackPts[k].size(); i=i+settings.samplingStep(0)){
	
	if ((trackPts[k][i].valid == 1) && (trackPts[k][i].DEMPt[0].valid == 1)){
	  Vector3 model;
	
	  model = trackPts[k][i].LOLAPt[2];
	  model.z() *= 1000; // LOLA data is in km, DEMGeo is in m (for LROC DTMs).
          
	  if ((model[0] >1e-100) && (model[1] > 1e-100) && (model[2]>1e-100) && 
              (model[0] < 360) && (model[0] > -180) && 
              (model[2] > 1720000) && (model[2] < 1750000)){//this should go into the LOLA reader 
	    if (verbose > 1){ 
	      cout<<"altitude="<<model[2]<<", dem="<< trackPts[k][i].DEMPt[0].val<<endl; 
	    }

	    Vector3 xyzModel = lon_lat_radius_to_xyz(model);
	    xyzModelArray.push_back(xyzModel);
	    llrModelArray.push_back(model);
            float radError = fabs(model[2]-1000*trackPts[k][i].DEMPt[0].val);
	    radErrorArray.push_back(radError);
            
	  }
	}
      }
    }
    
    xyzMatchArray.resize(xyzModelArray.size());
    xyzErrorArray.resize(xyzModelArray.size());
    
   
    if( verbose > 0 ){
      cout << "Number of points to be compared: " << xyzModelArray.size() << endl;
    }

    Vector3 modelCentroid;

    if ((settings.maxNumIter == 0) && (settings.matchWindowHalfSize(0)==1) && (settings.matchWindowHalfSize(1)==1)){
      //run error calculations with no re-estimation
      //GetMatchesFromDEM(xyzModelArray, llrModelArray, DEM, DEMGeo, settings.noDataVal, xyzMatchArray, xyzErrorArray);
    }
    else{
      //run ICP-matching
      modelCentroid = find_centroid(xyzModelArray);
      cout<<"modelCentroid="<<modelCentroid<<endl;
      ICP_LIDAR_2_DEM(xyzMatchArray, DEM, DEMGeo, xyzModelArray, llrModelArray, settings, 
		      currTranslation, currRotation, modelCentroid, xyzErrorArray);
    }
    string overlapDEMFileNoPathNoExt;
    cout<<"inputDEMFilename="<<inputDEMFilename<<endl;
    overlapDEMFileNoPathNoExt = GetFilenameNoExt(GetFilenameNoPath(inputDEMFilename));
    
    string LOLAFilenameNoPathNoExt;
    cout<<"inputCSVFilename="<<inputCSVFilename<<endl;
    LOLAFilenameNoPathNoExt = GetFilenameNoExt(GetFilenameNoPath(inputCSVFilename));
    cout<<"LOLAFilenameNoPathNoExt="<<LOLAFilenameNoPathNoExt<<endl;
    cout<<"overlapDEMFileNoPathNoExt="<<overlapDEMFileNoPathNoExt<<endl;
    string statsFilename;
    
    statsFilename = resDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapDEMFileNoPathNoExt + "_stats.txt";
    errorFilename = resDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapDEMFileNoPathNoExt + "_error.txt";
    
    if( vm.count("error-filename") && (xyzModelArray.size() >0 )){
      vector<string> titles(4);
      titles[0] = "Latitude";
      titles[1] = "Longitude";
      titles[2] = "Radius (m)";
      titles[3] = "Errors";  
      SaveDEMErrors( errorFilename, llrModelArray, xyzErrorArray, titles );
       
      vector<float> histBins;
      int numBins = 15;
      histBins.resize(numBins);
      histBins[0] = 0;
      histBins[1] = 25;
      histBins[2] = 50; 
      histBins[3] = 75;
      histBins[4] = 100;
      histBins[5] = 125;
      histBins[6] = 150; 
      histBins[7] = 175;
      histBins[8] = 200;
      histBins[9] = 225;
      histBins[10] = 250; 
      histBins[11] = 1000;
      histBins[12] = 2500;
      histBins[13] = 5000;
      histBins[14] = 10000; 
    
      SaveStatistics (statsFilename, /*xyzErrorArray*/radErrorArray, histBins);
    }
    
    if(( verbose >= 0 ) && (xyzModelArray.size() > 0)){ 

       /*
       //this must be changed - START
       ImageViewRef<PixelGray<float> >   interpDEM = interpolate(edge_extend(DEM.impl(),
       ConstantEdgeExtension()),
       BilinearInterpolation());
       //this must be changed - END
       */
   
      //int DEMcenterCol = interpDEM.cols() / 2;
      // int DEMcenterRow = interpDEM.rows() / 2;
      int DEMcenterCol = DEM.cols() / 2;
      int DEMcenterRow = DEM.rows() / 2;
      Vector2 lonlat = DEMGeo.pixel_to_lonlat( Vector2(DEMcenterCol,DEMcenterRow) );
      double DEMcenterR = DEMGeo.datum().radius(lonlat.x(),lonlat.y());
      //Vector3 DEMcenter_llr(lonlat.x(), lonlat.y(),interpDEM(DEMcenterCol,DEMcenterRow));
      Vector3 DEMcenter_llr(lonlat.x(), lonlat.y(),DEM(DEMcenterCol,DEMcenterRow));
      Vector3 DEMcenter_xyz = DEMGeo.datum().geodetic_to_cartesian(DEMcenter_llr);
      Vector3 translated_xyz = DEMcenter_xyz + currTranslation;
      Vector3 translated_llr = DEMGeo.datum().cartesian_to_geodetic(translated_xyz);
      Vector3 translation_llr = translated_llr - DEMcenter_llr;
      
      translation_llr.x() = DEMcenterR * tan( translation_llr.x() * M_PI/180.0  );
      translation_llr.y() = DEMcenterR * tan( translation_llr.y() * M_PI/180.0  );
      
      // Decompose rotation matrix to Euler Angles (phi, theta, psi about Z, X, and Z axes):
      Vector3 euler_angles = rotation_matrix_to_euler_xyz( currRotation );
      euler_angles *= 180/M_PI;
      
      // Decompose rotation matrix to Axis-angle:
      Vector3 axis_angle = matrix_to_axis_angle( currRotation );
      double axis_angle_deg = norm_2( axis_angle ) * 180/M_PI;
      
      // Work out estimate of the scale factor
      Vector3 f_cent = find_centroid( xyzMatchArray );
      
      Matrix<double> f_centered(xyzMatchArray.size(),3);
      Matrix<double> m_centered(xyzMatchArray.size(),3);
      
      for( unsigned int i = 0; i < xyzMatchArray.size(); i++ ){
	select_row(f_centered,i) = xyzMatchArray[i] - f_cent;
	select_row(m_centered,i) = xyzModelArray[i] - modelCentroid;
      }
      
      Matrix<double> ratio(xyzMatchArray.size(),3);
      Vector3 scale;
      
      ratio = elem_quot( f_centered * transpose(currRotation) , m_centered );
      for( unsigned int i = 0; i<=2; i++ ){
	scale(i) = sum( select_col(ratio,i) ) / ratio.rows();
      }
      
      cout << endl
	   << "Translation (xyz) = " << currTranslation << endl
	   << "Translation (llr) = " << translation_llr << endl
	   << "Rotation = " << currRotation << endl
	   << "Euler Angles (xyz in degrees)  = " << euler_angles << endl
	   << "Axis Angle = " << axis_angle << ", " << axis_angle_deg << " degrees" << endl
	   << "Scale factors = " << scale << endl
	//<< "Center = " << center << endl;
	   << "Centroid = " << modelCentroid << endl
	   << "allin1linehead: Translation xyz\tTraslation llr\tMatching Error\tEuler Angles xyz\tAxis Angle degrees\tScale Factors\tCentroid" << endl
	   << "alldatain1line: " << currTranslation << "\t" 
	   << translation_llr << "\t" 
	   << xyzErrorArray.sum()/xyzErrorArray.size() << "\t" 
	   << euler_angles << "\t" 
	   << axis_angle << " " << axis_angle_deg << "\t" 
	   << scale << "\t" 
	   << modelCentroid << endl;
    }
    
  }
  
  return 0;
}

