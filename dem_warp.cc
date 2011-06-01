#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography/GeoReference.h>
#include <vw/Cartography/PointImageManipulation.h>
#include <asp/Core/OrthoRasterizer.h>
using namespace vw;
using namespace vw::cartography;

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "DEMManipulation.h"

namespace vw {
  template<> struct PixelFormatID<Vector3> { static const PixelFormatEnum value=VW_PIXEL_GENERIC_3_CHANNEL; };
}

struct Options {
  // Input
  std::string input_name, matrix_string;

  // Settings
  float dem_spacing;
  double nodata_value;
  Matrix<double,3,4> affine;

  // Output
  std::string out_prefix;
};

class AffineFunc : public UnaryReturnSameType {
  Matrix<double,3,4> m_affine;
  Matrix<double,3,3> m_rotation;
  Vector3 m_translation;
public:
  AffineFunc( Matrix<double,3,4> const& affine ) : m_affine(affine) {
    m_rotation = submatrix(m_affine,0,0,3,3);
    m_translation = select_col(m_affine,3);
  }

  template <class T>
  T operator()( T const& p ) const {
    if ( p == T() ) return p;
    return m_rotation*p+m_translation;
  }
};

template <class ImageT>
UnaryPerPixelView<ImageT, AffineFunc>
inline rotate_pointcloud( ImageViewBase<ImageT> const& image,
                          Matrix<double,3,4> const& affine ) {
  return UnaryPerPixelView<ImageT,AffineFunc>( image.impl(),
                                               AffineFunc(affine) );
}

void handle_arguments( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("");
  general_options.add_options()
    ("dem-spacing,s", po::value(&opt.dem_spacing)->default_value(0.0), "Set the DEM post size (if this value is 0, the post spacing size is computed for you)")
    ("output-prefix,o", po::value(&opt.out_prefix), "Specify the output prefix")
    ("nodata-value", po::value(&opt.nodata_value)->default_value(std::numeric_limits<double>::max() ), "Nodata value used in the input DEM and also in the output DEM. Will auto-detect if possible.")
    ("help,h", "Display this help message");

  po::options_description positional("");
  positional.add_options()
    ("input-sub", po::value(&opt.input_name), "Subject DEM")
    ("input-matrix", po::value(&opt.matrix_string), "Input Matrix");

  po::positional_options_description positional_desc;
  positional_desc.add("input-sub",    1);
  positional_desc.add("input-matrix", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " [options] <input dem> <matrix string> \n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.input_name.empty() ||
       opt.matrix_string.empty() )
    vw_throw( ArgumentErr() << "Missing one of the input arguments!\n"
              << usage.str() << general_options );
}

int main( int argc, char *argv[] ) {
  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    cartography::GeoReference georef;
    if ( !read_georeference( georef, opt.input_name ) )
      vw_throw( ArgumentErr() << "Failed to read input file's georeference." );
    DiskImageView<double> dem( opt.input_name );

    { // Attempt to extract nodata value
      DiskImageResource *disk_dem_rsrc =
        DiskImageResource::open(opt.input_name);
      if ( opt.nodata_value != std::numeric_limits<double>::max() ) {
        vw_out() << "\t--> Using user-supplied nodata value: " << opt.nodata_value << ".\n";
      } else if ( disk_dem_rsrc->has_nodata_read() ) {
        opt.nodata_value = disk_dem_rsrc->nodata_read();
        vw_out() << "\t--> Extracted nodata value from file: " << opt.nodata_value << ".\n";
      } else {
        opt.nodata_value = min_pixel_value( dem );
        vw_out() << "\t--> Using min value of DEM for nodata value: " << opt.nodata_value << ".\n";
      }
    }

    // Parse out the matrix transform
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(",:;{}[] ");
    tokenizer tokens( opt.matrix_string, sep );
    Matrix<double,3,4>::iterator mat_iter = opt.affine.begin();
    for ( tokenizer::iterator tok_iter = tokens.begin();
          tok_iter != tokens.end(); tok_iter++ ) {
      if ( mat_iter == opt.affine.end() )
        vw_throw( ArgumentErr() << "Input matrix string has incorrect number of elements. I need a 3x4 matrix." );
      *mat_iter = boost::lexical_cast<double>(*tok_iter);
      mat_iter++;
    }
    if ( mat_iter != opt.affine.end() )
      vw_throw( ArgumentErr() << "Input matrix string has incorrect number of elements. I need a 3x4 matrix." );
    vw_out() << "\t--> Matrix: " << opt.affine << "\n";

    // Convert to rotated Point Cloud
    ImageViewRef<Vector3> point_ref =
      rotate_pointcloud( dem_to_point_cloud(create_mask(dem, opt.nodata_value),
                                            georef), opt.affine );
    DiskCacheImageView<Vector3> point_cache( point_ref, "tif",
                                             TerminalProgressCallback("dem","Cache: ") );

    // Convert back to a projected image after applying distortion
    int32 subsample_amt = int32(norm_2(Vector2i(point_ref.cols(),point_ref.rows()))) / 1024;
    if (subsample_amt < 1 ) subsample_amt = 1;
    Vector3 avg_location =
      mean_pixel_value(subsample(point_cache, subsample_amt));
    ImageViewRef<Vector3> projected_ref =
      project_point_image(xyz_to_lon_lat_radius(point_cache,true,avg_location[0] >= 0 ), georef );

    // Rasterize back into a DEM
    OrthoRasterizerView<PixelGray<float> >
      rasterizer(projected_ref,
                 select_channel(projected_ref,2),
                 opt.dem_spacing);
    rasterizer.set_use_minz_as_default(false);
    rasterizer.set_default_value(opt.nodata_value);

    BBox3 dem_bbox = rasterizer.bounding_box();
    vw_out() << "\nDEM Bounding box: " << dem_bbox << "\n";

    Matrix3x3 georef_affine_transform = rasterizer.geo_transform();
    vw_out() << "Georeferencing Transform: " << georef_affine_transform << "\n";
    cartography::GeoReference georef_out = georef;
    georef_out.set_transform(georef_affine_transform);

    // Work out what the output prefix should be
    if ( opt.out_prefix.empty() ) {
      int index = opt.input_name.rfind(".");
      opt.out_prefix = opt.input_name;
      if ( index != -1 )
        opt.out_prefix.erase(index, opt.out_prefix.size());
    }

    vw_out() << "\nWriting DEM.\n";
    ImageViewRef<PixelGray<float> > block_dem_raster =
      block_cache(rasterizer, Vector2i(rasterizer.cols(), 128), 0);
    DiskImageResourceGDAL rsrc( opt.out_prefix + "-DEM.tif",
                                block_dem_raster.format() );
    rsrc.set_nodata_write( opt.nodata_value );
    write_georeference( rsrc, georef );
    write_image( rsrc, block_dem_raster, TerminalProgressCallback("dem","") );

  } catch ( const ArgumentErr& e ) {
    vw_out() << e.what() << std::endl;
    return 1;
  } catch ( const Exception& e ) {
    vw_out() << "VW Error: " << e.what() << std::endl;
    return 1;
  } catch ( std::bad_alloc& e ) {
    std::cerr << "Error: Ran out of Memory!" << std::endl;
    return 1;
  }
}
