#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// Erases a file path if one exists and returns the base string 
std::string GetFilenameNoPath(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind("/");
  if (index != -1)
    result.erase(0, index+1);
  return result;
}

/// Erases a file suffix if one exists and returns the base string
std::string GetFilenameNoExt(std::string const& filename) {
  std::string result = filename;
  int index = result.rfind(".");
  if (index != -1)
      result.erase(index, result.size());
  return result;
}

static bool readStringListFromFile(const std::string& filename,
                                   std::vector<std::string>& l)
{
  std::ifstream file;
  file.open(filename.c_str());
  char str[255];
  while ( !file.eof() ) {
     std::string temp;
     file.getline(str, 255);
     l.push_back(std::string(str));
     std::cout<<str<<std::endl;
 }

  return true;
file.close();

}
/*
static bool readStringList( const std::string& filename,
                            std::vector<std::string>& l ) {
  l.resize(0);
  cv::FileStorage fs(filename, cv::FileStorage::READ);
  if( !fs.isOpened() )
    return false;
  cv::FileNode n = fs.getFirstTopLevelNode();
  if( n.type() != cv::FileNode::SEQ )
    return false;
  cv::FileNodeIterator it = n.begin(), it_end = n.end();
  for( ; it != it_end; ++it )
    l.push_back((std::string)*it);
  return true;
}
*/
void convertToJpg( std::vector<std::string> image_list)
{
  for (int i=0; i < image_list.size(); i++){
    std::cout<<image_list[i]<<std::endl;
    cv::Mat image = cv::imread(image_list[i], 0);
    std::string newFilename = GetFilenameNoExt(image_list[i])+".jpg";
     std::cout<<newFilename<<std::endl;   
     cv::imwrite(newFilename, image);
  }
} 

static void calculate_chessboard_corners( cv::Size const& board_size, double square_size, std::vector<cv::Point3f>& corners ) {
  corners.resize(0);
  for( int i = 0; i < board_size.height; i++ )
    for( int j = 0; j < board_size.width; j++ )
      corners.push_back(cv::Point3f(float(j*square_size),
                                    float(i*square_size), 0));
}

std::ostream& clarity_print( std::ostream& out, cv::Mat const& mat ) {
  out << "[ ";
  for ( int row = 0; row < mat.rows; row++ ) {
    for ( int col = 0; col < mat.cols; col++ ) {
      if ( row != 0 || col != 0 )
        out << " ";
      out << mat.at<double>(row,col);
    }
  }
  out << "]";
  return out;
}

int main( int argc, char* argv[] ) {
  cv::Size board_size, image_size;
  std::string input_file;
  double square_size;

  po::options_description general_options("Options");
  general_options.add_options()
    ("width,w", po::value(&board_size.width), "Number of inside corners, width-wise")
    ("height,h", po::value(&board_size.height), "Number of inside corners, height-wise")
    ("square_size,s", po::value(&square_size)->default_value(1), "Square size in whatever units you want to use.")
    ("show-corners", "Show corners detected automatically")
    ("show-rectified", "Show rectified images")
    ("help", "Display this help message");

  po::options_description hidden_options("");
  hidden_options.add_options()
    ("input-file", po::value(&input_file));

  po::options_description options("Allowed Options");
  options.add(general_options).add(hidden_options);

  po::positional_options_description p;
  p.add("input-file", 1);

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input yaml>\n\n";
  usage << general_options << std::endl;

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(options).positional(p).run(), vm );
    po::notify( vm );
  } catch ( const po::error& e ) {
    std::cout << "Error parsing command line arguments.\n";
    std::cout << "\t" << e.what() << "\n\n";
    std::cout << usage.str();
    return 1;
  }
  if ( vm.count("help" ) ) {
    std::cout << usage.str();
    return 1;
  }

  std::vector<std::string> image_list;
  //if ( !readStringList(input_file,image_list) ) {
  if ( !readStringListFromFile(input_file,image_list) ) {
    std::cout << "Unable to open \"" << input_file << "\".\n";
    std::cout << usage.str();
    return 1;
  }

  // 0.) convert to JPG
  //convertToJpg(image_list);
 
  // 1.) Loading up images and finding checkobard pattern
  std::cout << "Loading images:\n";
  std::vector<std::vector<cv::Point2f> > image_measurements[2];
  image_measurements[0].resize(image_list.size()/2);
  image_measurements[1].resize(image_list.size()/2);
  std::cout<<"num image pairs = "<< image_list.size()/2<<std::endl;
  for ( size_t pair_idx = 0; pair_idx < image_list.size()/2; pair_idx++ ) {
    bool pair_good = true;
    for ( size_t j = 0; j < 2 && pair_good; j++ ) {
      cv::Mat image = cv::imread( image_list[2*pair_idx+j], 0 );
      if ( pair_idx == 0 && j == 0 )
        image_size = image.size();
      pair_good &=
        cv::findChessboardCorners( image, board_size,
                                   image_measurements[j][pair_idx],
                                   CV_CALIB_CB_ADAPTIVE_THRESH +
                                   CV_CALIB_CB_NORMALIZE_IMAGE );
      if ( !pair_good ) {
        std::cout << "Rejected: " << image_list[2*pair_idx+j] << std::endl;
        break;
      }
      cv::cornerSubPix( image, image_measurements[j][pair_idx],
                        cv::Size(11,11), cv::Size(-1,-1),
                        cv::TermCriteria( CV_TERMCRIT_EPS + CV_TERMCRIT_ITER,
                                          30, 0.1 ));
      if ( vm.count("show-corners") ) {
        std::cout << image_list[2*pair_idx+j] << "\n";
        cv::Mat cimg;
        cvtColor(image, cimg, CV_GRAY2BGR);
        cv::drawChessboardCorners(cimg, board_size,
                                  image_measurements[j][pair_idx], true);
	
        imshow("corners", cimg);
        cv::waitKey(1);
        //char c;
        //std::cin >> c;
      } else {
        std::cout << "." << std::flush;
      }
    } // end j loop
    if ( !pair_good ) {
      image_measurements[0][pair_idx].resize(0);
      image_measurements[1][pair_idx].resize(0);
    }
  }
  std::cout << "\n\n";

  // 1.a) Assign positions for our captured corners
  std::vector<std::vector<cv::Point3f> > object_points(1);
  calculate_chessboard_corners(board_size, square_size, object_points[0]);
  object_points.resize(image_measurements[0].size(),object_points[0]);

  // 2.) Solve for individual camera's intrinsic parameters.
  // #if 0
  cv::Mat camera_matrix[2], dist_coeffs[2];
  for ( size_t j = 0; j < 2; j++ ) {
    camera_matrix[j] = cv::Mat::eye(3, 3, CV_64F);
    dist_coeffs[j] = cv::Mat::zeros(8, 1, CV_64F);
    //#if 0
    std::vector<cv::Mat> rvecs, tvecs;
     
    double rms;  

     std::cout << "init camera[" << j << "]" << "\n";
     //std::cout<<"dist_coeffs"<<dist_coeffs[j]<<"\n";
     //#if 0
    rms =
      cv::calibrateCamera(object_points, image_measurements[j], image_size,
                          camera_matrix[j], dist_coeffs[j], rvecs, tvecs//,
                          //CV_CALIB_USE_INTRINSIC_GUESS |
                          //CV_CALIB_ZERO_TANGENT_DIST | 
			  /*
                          CV_CALIB_FIX_K3 | CV_CALIB_FIX_K4 | 
                          CV_CALIB_FIX_K5 | CV_CALIB_FIX_K6*/ );
    
    std::cout << "camera[" << j << "] RMS error 2nd pass: " << rms << "\n";
    //std::cout<<"rvecs["<<j<<"]="<<rvecs[j]<<" "<<tvecs[j]<<endl;
  
    //dist_coeffs[j] = dist_coeffs[j](cv::Range(0,5), cv::Range::all() );

    std::cout << "camera_matrix_step1[" <<j<< "]=" << camera_matrix[j] << "\n"
	      << "  dist_coeffs_step1[" <<j<< "]=" <<dist_coeffs[j] << "\n";
    //#endif
  }

  // 3.) Solve for their extrinsic line up.
  cv::Mat R, T, E, F;
  //cv::Mat camera_matrix[2], dist_coeffs[2];
  std::vector<cv::Mat> rvecs, tvecs;
  double rms =
    cv::stereoCalibrate( object_points, image_measurements[0],
                         image_measurements[1],
                         camera_matrix[0], dist_coeffs[0],
                         camera_matrix[1], dist_coeffs[1],
                         image_size, R, T, E, F,
                         cv::TermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, 30, 1e-5 ),
			 /*		 
                         CV_CALIB_FIX_ASPECT_RATIO +
			 CV_CALIB_ZERO_TANGENT_DIST +
			 CV_CALIB_SAME_FOCAL_LENGTH +
			 CV_CALIB_RATIONAL_MODEL +
			 CV_CALIB_FIX_K3 + CV_CALIB_FIX_K4 + CV_CALIB_FIX_K5);
			 */
  
                         
                         CV_CALIB_USE_INTRINSIC_GUESS + //Ara 04/23
			 //CV_CALIB_FIX_ASPECT_RATIO + //Ara 04/23
			 //CV_CALIB_FIX_FOCAL_LENGTH + //Ara 04/23
			 //CV_CALIB_ZERO_TANGENT_DIST+
			 //CV_CALIB_SAME_FOCAL_LENGTH + //Ara 04/23
			 CV_CALIB_RATIONAL_MODEL +
			 CV_CALIB_FIX_K3 + CV_CALIB_FIX_K4 + CV_CALIB_FIX_K5);
 		 
  
  //cout << "done with RMS error=" << rms << endl;

  std::cout << "stereo calibrate RMS error: " << rms << "\n";

  dist_coeffs[0] = dist_coeffs[0](cv::Range(0,4), cv::Range::all() );
  dist_coeffs[1] = dist_coeffs[1](cv::Range(0,4), cv::Range::all() );

  std::cout << "camera_matrix[0] " << camera_matrix[0] << "\n"
            << "  dist_coeffs[0] " << dist_coeffs[0] << "\n"
            << "camera_matrix[1] " << camera_matrix[1] << "\n"
            << "  dist_coeffs[1] " << dist_coeffs[1] << "\n\n";

  std::cout << "R: " << R << "\n";
  std::cout << "T: " << T << "\n";

  //3.bis) 
  // CALIBRATION QUALITY CHECK
  // because the output fundamental matrix implicitly
  // includes all the output information,
  // we can check the quality of calibration using the	/
  // epipolar geometry constraint: m2^t*F*m1=0
  double err = 0;
  int npoints = 0;
  /*
  vector<Vec3f> lines[2];
  
 for( i = 0; i < nimages; i++ ){
    
    int npt = (int)imagePoints[0][i].size();
    Mat imgpt[2];
    for( k = 0; k < 2; k++ ){
      imgpt[k] = Mat(imagePoints[k][i]);
      undistortPoints(imgpt[k], imgpt[k], cameraMatrix[k], distCoeffs[k], Mat(), cameraMatrix[k]);
      computeCorrespondEpilines(imgpt[k], k+1, F, lines[k]);
    }
    for( j = 0; j < npt; j++ ){
      double errij = fabs(imagePoints[0][i][j].x*lines[1][j][0] +
			  imagePoints[0][i][j].y*lines[1][j][1] + lines[1][j][2]) +
	fabs(imagePoints[1][i][j].x*lines[0][j][0] +
	     imagePoints[1][i][j].y*lines[0][j][1] + lines[0][j][2]);
      err += errij;
    }
    npoints += npt;
  }
  cout << "average reprojection err = " <<  err/npoints << endl;
  */

  // 4.) Write out the solution in Clarity format
  {
    std::ofstream f("irg_stereo_calib.cameras");
    f << std::setprecision(9) << "\n";
    // 4.1) Write left
    f << "{ type=Camera_Model_Matrix; camera_model_version=1; name=stereo_left;\n"
      << "  height=" << image_size.height << "; width=" << image_size.width << ";\n"
      << "  frame={ ver=1; name=stereoLeft; parent=linkageFrame; wrt=;\n"
      << "          loc={ x=0; y=0; z=0; rx=0; ry=0; rz=0; }; };\n"
      << "  version=2;\n"
      << "  distortion={ version=1;\n"
      << "               elements=";
    clarity_print( f, dist_coeffs[0] ) << ";\n"
      << "               orientation=r; };\n"
      << "  matrix={ nrows=3; ncols=3;\n"
      << "           elements=";
    clarity_print( f, camera_matrix[0] ) << "; };\n}\n";

    // 4.2) Solve for the right to left frame store parameters
    cv::Mat R_prime = R.inv();
    cv::Mat T_prime = -R_prime*T;
    double pitch = std::atan2(-R_prime.at<double>(2,0),
                              sqrt(R_prime.at<double>(0,0)*R_prime.at<double>(0,0) +
                                   R_prime.at<double>(1,0)*R_prime.at<double>(1,0)));
    double roll = 0, yaw = 0;
    // This code is copied from clarity's rotation_matrix.h just to be safe
    if ( fabs( pitch-0.5*M_PI ) < std::numeric_limits<double>::epsilon() ) {
      roll = std::atan2(R_prime.at<double>(0,1), R_prime.at<double>(1,1));
    }
    if ( fabs( pitch+0.5*M_PI ) < std::numeric_limits<double>::epsilon() ) {
      roll = -std::atan2(R_prime.at<double>(0,1), R_prime.at<double>(1,1));
    } else {
      roll = std::atan2(R_prime.at<double>(2,1), R_prime.at<double>(2,2));
      yaw  = std::atan2(R_prime.at<double>(1,0), R_prime.at<double>(0,0));
    }

    // 4.3) Write right
    f << "{ type=Camera_Model_Matrix; camera_model_version=1; name=stereo_right;\n"
      << "  height=" << image_size.height << "; width=" << image_size.width << ";\n"
      << "  frame={ ver=1; name=stereoRight; parent=stereoLeft; wrt=;\n"
      << "          loc={ x=" << T_prime.at<double>(0,0)
                    << "; y=" << T_prime.at<double>(0,1)
                    << "; z=" << T_prime.at<double>(0,2)
      << "; rx=" << roll << "; ry=" << pitch << "; rz=" << yaw << "; }; };\n"
      << "  version=2;\n"
      << "  distortion={ version=1;\n"
      << "               elements=";
    clarity_print( f, dist_coeffs[1] ) << ";\n"
      << "               orientation=r; };\n"
      << "  matrix={ nrows=3; ncols=3;\n"
      << "           elements=";
    clarity_print( f, camera_matrix[1] ) << "; };\n}\n";
  }

  //4.bis save in txt format (non clarity)
  std::ofstream f("stereo_calibration.txt");
  f<<"CAMERA_MATRIX_LEFT: "<<camera_matrix[0].at<double>(0,0)<<" "<<camera_matrix[0].at<double>(0,1)<<" "<<camera_matrix[0].at<double>(0,2)
                           <<" "<<camera_matrix[0].at<double>(1,0)<<" "<<camera_matrix[0].at<double>(1,1)<<" "<<camera_matrix[0].at<double>(1,2)
                           <<" "<<camera_matrix[0].at<double>(2,0)<<" "<<camera_matrix[0].at<double>(2,1)<<" "<<camera_matrix[0].at<double>(2,2)<<std::endl;
  f<<"DISTORSION_COEFFICIENTS_LEFT: "<<dist_coeffs[0].at<double>(0,0)<<" "<<dist_coeffs[0].at<double>(0,1)<<" "<<dist_coeffs[0].at<double>(0,2)<<" "<<dist_coeffs[0].at<double>(0,4)<<std::endl;
  f<<"CAMERA_MATRIX_RIGHT: "<<camera_matrix[1].at<double>(0,0)<<" "<<camera_matrix[1].at<double>(0,1)<<" "<<camera_matrix[1].at<double>(0,2)
                           <<" "<<camera_matrix[1].at<double>(1,0)<<" "<<camera_matrix[1].at<double>(1,1)<<" "<<camera_matrix[1].at<double>(1,2)
                           <<" "<<camera_matrix[1].at<double>(2,0)<<" "<<camera_matrix[1].at<double>(2,1)<<" "<<camera_matrix[1].at<double>(2,2)<<std::endl;
  f<<"DISTORSION_COEFFICIENTS_RIGHT: "<<dist_coeffs[1].at<double>(0,0)<<" "<<dist_coeffs[1].at<double>(0,1)<<" "<<dist_coeffs[1].at<double>(0,2)<<" "<<dist_coeffs[1].at<double>(0,4)<<std::endl;
  f<<"ROTATION_MATRIX: "<<R.at<double>(0,0)<<" "<<R.at<double>(0,1)<<" "<<R.at<double>(0,2)
                        <<" "<<R.at<double>(1,0)<<" "<<R.at<double>(1,1)<<" "<<R.at<double>(1,2)
                        <<" "<<R.at<double>(2,0)<<" "<<R.at<double>(2,1)<<" "<<R.at<double>(2,2)<<std::endl;
  f<<"TRANSLATION_VECTOR: "<<T.at<double>(0,0)<<" "<<T.at<double>(0,1)<<" "<<T.at<double>(0,2)<<std::endl;
  f<<"IMAGE_WIDTH_HEIGHT: "<<image_size.width<<" "<<image_size.height<<std::endl;
  f.close();

  
  // 5.) Show off the solution with epirectified images
  if ( vm.count("show-rectified") ) {
    // 5.a) calculate ideal rotation and new P's for epipolar
    // rectification.
    cv::Mat R1, R2, P1, P2, Q;
    cv::Rect validRoi[2];
  
    cv::stereoRectify( camera_matrix[0], dist_coeffs[0],
                       camera_matrix[1], dist_coeffs[1],
                       image_size, R, T, 
                       R1, R2, P1, P2, 
                       Q,  CV_CALIB_ZERO_DISPARITY, 1,
                       image_size, &validRoi[0], &validRoi[1]);
               
    
    // 5.b) calculate the rectification maps
    cv::Mat rmap[2][2];
    cv::initUndistortRectifyMap(camera_matrix[0], dist_coeffs[0], R1, P1, image_size, CV_16SC2, rmap[0][0], rmap[0][1]);
    cv::initUndistortRectifyMap(camera_matrix[1], dist_coeffs[1], R2, P2, image_size, CV_16SC2, rmap[1][0], rmap[1][1]);
    
    // 5.c) rendering
    cv::Mat canvas;
    cv::Size canvas_sub_size;
    double sf = 1;//300.0 / std::max(image_size.width,image_size.height);
    canvas_sub_size.width = image_size.width;// * sf + 0.5;
    canvas_sub_size.height = image_size.height;// * sf + 0.5;
    canvas.create(canvas_sub_size.height,canvas_sub_size.width*2,CV_8UC3);
    for ( size_t i = 0; i < image_measurements[0].size(); i++ ) {
      for ( size_t j = 0; j < 2; j++ ) {
        cv::Mat img = cv::imread(image_list[i*2+j], 0), rimg, cimg;
        cv::remap(img, rimg, rmap[j][0], rmap[j][1], CV_INTER_LINEAR);
        cv::cvtColor(rimg, cimg, CV_GRAY2BGR);
        cv::Mat canvas_part =  // Lazy copy?
          canvas(cv::Rect(canvas_sub_size.width*j, 0,
                          canvas_sub_size.width, canvas_sub_size.height));
        cv::resize(cimg, canvas_part, canvas_part.size(), 0, 0, CV_INTER_AREA);
      }
      // Drawing lines to help identify correctness
      for ( size_t j = 0; j < canvas.rows; j += 16 )
        line( canvas, cv::Point(0,j), cv::Point(canvas.cols, j), cv::Scalar(0,255,0), 1, 8 );
      cv::imshow("rectified", canvas);
      cv::waitKey(500);
    }
    
  }
  
  return 0;
}
