// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

// #include <iostream>
// #include <string>
// #include <fstream>
// #include <vector>
// #include <sys/types.h>
// #include <sys/stat.h>
// #include <math.h>

//boost
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <boost/operators.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
namespace fs = boost::filesystem;

//VW
// #include <vw/Core.h>
// #include <vw/Image.h>
// #include <vw/FileIO.h>
// #include <vw/Cartography.h>
#include <vw/Math/EulerAngles.h>

//ATK
#include "coregister.h"
#include "icp.h"
#include "tracks.h"
#include "util.h"

using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

int main( int argc, char *argv[] )
{
 int verbose;
 std::vector<std::string> inputDEMFiles;
 std::string precDEMDir="NO_DIR";

 po::options_description general_options("Options");
 general_options.add_options()
    ("lidar-filename,l", po::value<std::string>())
    ("inputDEMFiles,d", po::value<std::vector<std::string> >(&inputDEMFiles))
    ("DEM-precision-directory,p", 
     po::value<std::string>(&precDEMDir), 
		"DEM precision directory.") 
    ("results-directory,r", 
		po::value<std::string>()->default_value("./results"), 
		"results directory.") 
    ("settings-filename,s", 
		po::value<std::string>()->default_value("lidar2dem_settings.txt"), 
		"settings filename.")
    ("errors,e", 
		po::value<bool>()->zero_tokens(), 
		"Write out detailed errors.")
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
  try{
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  }
  catch( po::error& e ){
    std::cout << "An error occured while parsing command line arguments." << endl;
    std::cout << "\t" << e.what() << endl << endl;
    std::cout << usage.str();
    return 1;
  }
  
  if( vm.count("help") )
    {
      std::cerr << usage.str() << std::endl;
      return 0;
    }

  if( !vm.count("lidar-filename") ){
      std::cerr << "Error: Must specify at least  one Lidar file!" << std::endl << std::endl;
      std::cerr << usage.str();
      return 1;
  }

  if( inputDEMFiles.empty() ) {
    std::cerr << "Error: Must specify at least  one DEM file!" << std::endl << std::endl;
    std::cerr << usage.str();
    return 1;
  }

  //determine if inputDEMFiles is a text file containing a list of DEM files, one DEM file or a set of DEM files
  //by reading the file extension. A text file containing the DEM list *must* have extension .txt 
  vector<string> DEMFiles( AccessDataFilesFromInput(inputDEMFiles) );
 
  //Set up VW logging
  vw_log().console_log().rule_set().add_rule(vw::InfoMessage,"*");

 
  struct CoregistrationParams settings;
  if( ReadConfigFile( vm["settings-filename"].as<string>(), &settings) && verbose > 0 ){
   cout << "Config file " << vm["settings-filename"].as<string>() << " found." << endl;
  } 
  else if( verbose > 0 ){
   cout << "Config file " << vm["settings-filename"].as<string>() << " not found, using defaults." << endl;
  }

  if( verbose > 0 ){ 
    cout << settings << endl;
    if( DEMFiles.size() > 1 ){
      string joined = boost::algorithm::join(DEMFiles, ",");
      cout << DEMFiles.size() << " DEM files( " << joined << " )" << endl;
    }
    else{ cout << "DEM file: " << DEMFiles[0] << endl; }
  }

  fs::path results_path( vm["results-directory"].as<string>() );
  fs::path aux_path( results_path/"aux" );

  if( verbose > 0 ){ cout << "aux_path= " << aux_path << endl; }
  //create the results and auxiliary directories
  if( fs::create_directory( results_path ) == false ){
    cerr << "Could not create the " << results_path << " directory." << endl;
    return 1;
  }
  if( fs::create_directory( aux_path ) == false ){
    cerr << "Could not create the " << aux_path << " directory." << endl;
    return 1;
  }


  //read LOLA tracks
  vector<vector<LOLAShot> > trackPts;
  fs::path tracks_path( vm["lidar-filename"].as<string>() );
  try{
    trackPts = LOLAFileRead( tracks_path.string() );
    if( verbose > 0 ){ cout << "Tracks file: " << tracks_path.string() << endl; }
  }
  catch( const vw::IOErr& error ){
    cerr << error.what() << endl;
    exit(1);
  }
 
  //determine the overlapping DEMs - START
  cout << "Selecting the overlapping DEMs ..." << endl;
  std::vector<int> overlapIndices;
  fs::path overlapList_path( aux_path/tracks_path.filename().replace_extension("_overlap_list.txt") );
  try{ 
    overlapIndices = ReadVectorFrom<int>(overlapList_path);
    if( !overlapIndices.empty() and overlapIndices[0] == -1 ){
      // The previous generation of this code could have written a -1 to the file.
      overlapIndices.clear();
      }
  }
  catch( const vw::ArgumentErr& error){
    //Vector4 lat_lon_bb = FindMinMaxLat(trackPts); 
    BBox2 bb = FindShotBounds( trackPts );
    overlapIndices = makeOverlapList( DEMFiles, bb );
    if( !overlapIndices.empty() ){
      WriteVectorTo( overlapList_path, overlapIndices );
    }
  }
  if( verbose > 0 ){
    PrintOverlapList( overlapIndices );
    cout << "done determining the overlapping DEMs." << endl;
  }
  //determine the overlapping DEMs - END
 
  for( unsigned int index = 0; index < overlapIndices.size(); ++index ){
    //string inputDEMFilename = DEMFiles[overlapIndices[index]];
    fs::path inputDEM_path( DEMFiles[overlapIndices[index]] );

    if( verbose > 0 ){ 
      cout << "DEM filename: " << inputDEM_path << endl; 
    }
    
    //read DEM file      
    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceGDAL( inputDEM_path.string() ) );
    if( rsrc->has_nodata_read() ){
      settings.noDataVal = rsrc->nodata_read();
	  if( verbose > 0 ){
        cout << "Found nodata value, " << settings.noDataVal 
             << ", in " << inputDEM_path << endl;
      }
    }
    else if( verbose > 0 ){
      cout << "Using default nodata value: " << settings.noDataVal << endl;
    }
    
    DiskImageView<float> DEM(rsrc);
    GeoReference DEMGeo;
    read_georeference(DEMGeo, inputDEM_path.string() );

    if( precDEMDir.compare(string("NO_DIR"))!=0 ){
      //select DEM points with high precision closest to LOLA tracks
      fs::path precdir_path( precDEMDir );
      fs::path prec_path( precdir_path/
                          boost::algorithm::replace_all_copy(inputDEM_path.string(), "dem", "Prec") );
      if( verbose > 0 ){ cout << "PrecFilename=" << prec_path << endl; }

      boost::shared_ptr<DiskImageResource> prec_rsrc( new DiskImageResourceGDAL( prec_path.string() ) );
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
    currRotation.set_identity();
    vector<Vector3> xyzModelArray;//LOLA
    vector<Vector3> llrModelArray;//LOLA
    vector<Vector3> xyzMatchArray;//DEM
    valarray<float> xyzErrorArray;//error array same size as model and feature
    vector<float>   radErrorArray; //altitude-radial error array
    
    //copy info to featureArray and modelArray
    for( unsigned int k = 0; k < trackPts.size(); ++k ){
      for( unsigned int i = 0; i < trackPts[k].size(); i += settings.samplingStep(0) ){
        if( (trackPts[k][i].valid == 1) && 
            (trackPts[k][i].DEMPt[0].valid) &&
            (trackPts[k][i].LOLAPt.size() > 2) ){
          Vector3 model( trackPts[k][i].LOLAPt[2] );
          model.z() *= 1000; // LOLA data is in km, DEMGeo is in m (for LROC DTMs).
          
          if( verbose > 1 ){ 
            cout << "altitude=" << model.z() << ", dem=" << trackPts[k][i].DEMPt[0].val << endl; 
          }
          Vector3 xyzModel = lon_lat_radius_to_xyz(model);
          xyzModelArray.push_back(xyzModel);
          llrModelArray.push_back(model);
          float radError = fabs(model[2]-1000*trackPts[k][i].DEMPt[0].val);
          radErrorArray.push_back(radError);
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
      cout << "modelCentroid=" << modelCentroid << endl;
      ICP_LIDAR_2_DEM( xyzMatchArray, DEM, DEMGeo, xyzModelArray, llrModelArray, settings, 
                       currTranslation, currRotation, modelCentroid, xyzErrorArray);
    }

    //string overlapDEMFileNoPathNoExt;
    //cout<<"inputDEMFilename="<<inputDEMFilename<<endl;
    //overlapDEMFileNoPathNoExt = GetFilenameNoExt(GetFilenameNoPath(inputDEMFilename));
    //overlapDEMFileNoPathNoExt = inputDEM_path.filename().stem();
    
    //string LOLAFilenameNoPathNoExt;
    //cout<<"inputCSVFilename="<<inputCSVFilename<<endl;
    //LOLAFilenameNoPathNoExt = GetFilenameNoExt(GetFilenameNoPath(inputCSVFilename));
    //LOLAFilenameNoPathNoExt = tracks_path.filename().stem();
    //cout<<"LOLAFilenameNoPathNoExt="<<LOLAFilenameNoPathNoExt<<endl;
    //cout<<"overlapDEMFileNoPathNoExt="<<overlapDEMFileNoPathNoExt<<endl;
    //string statsFilename;
    
    //statsFilename = resDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapDEMFileNoPathNoExt + "_stats.txt";
    ostringstream os;
    os << tracks_path.filename().stem().string() << "_" << inputDEM_path.filename().stem().string();
    fs::path stats_path( results_path/os.str().append("_stats.txt") );

    //errorFilename = resDir + "/" + LOLAFilenameNoPathNoExt + "_" + overlapDEMFileNoPathNoExt + "_error.txt";
    if( vm.count("errors") && !xyzModelArray.empty() ){
      fs::path error_path( results_path/os.str().append("_error.txt") );
      vector<string> titles(4);
      titles[0] = "Latitude";
      titles[1] = "Longitude";
      titles[2] = "Radius (m)";
      titles[3] = "Errors";  
      SaveDEMErrors( error_path.string(), llrModelArray, xyzErrorArray, titles );
       
      vector<float> histBins( 15 );
      for( unsigned int i = 0; i < 11; ++i ){
        histBins[i] = i * 25;
      }
      histBins[11] = 1000;
      histBins[12] = 2500;
      histBins[13] = 5000;
      histBins[14] = 10000; 
    
      SaveStatistics( stats_path.string(), radErrorArray, histBins );
    }
    
    if( ( verbose >= 0 ) && !xyzModelArray.empty() ){ 

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
      
      for( unsigned int i = 0; i < xyzMatchArray.size(); ++i ){
        select_row(f_centered,i) = xyzMatchArray[i] - f_cent;
        select_row(m_centered,i) = xyzModelArray[i] - modelCentroid;
      }
      
      Matrix<double> ratio(xyzMatchArray.size(),3);
      Vector3 scale;
      
      ratio = elem_quot( f_centered * transpose(currRotation) , m_centered );
      for( unsigned int i = 0; i<=2; ++i ){
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

