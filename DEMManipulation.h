#ifndef __DEMA_DEM_MANIPULATION_H__
#define __DEMA_DEM_MANIPULATION_H__

#include <vw/Cartography/PointImageManipulation.h>
#include <vw/Cartography/Datum.h>

namespace vw {

  class LLAtoXYZFunctor : public ReturnFixedType<Vector3> {
    cartography::Datum m_datum;
  public:
    LLAtoXYZFunctor( cartography::Datum const& datum ) : m_datum(datum) {}

    Vector3 operator()( Vector3 const& lla ) const {
      if ( lla == Vector3() )
        return Vector3();
      return m_datum.geodetic_to_cartesian( lla );
    }
  };

  template <class ImageT>
  UnaryPerPixelView<BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, cartography::DemToPointImageFunctor<typename ImageT::pixel_type> >, LLAtoXYZFunctor>
  dem_to_point_cloud( ImageViewBase<ImageT> const& image,
                      cartography::GeoReference const& georef ) {
    typedef cartography::DemToPointImageFunctor<typename ImageT::pixel_type> func1_type;
    typedef LLAtoXYZFunctor func2_type;
    typedef BinaryPerPixelView<PerPixelIndexView<VectorIndexFunctor>, ImageT, func1_type> inner_view;
    return UnaryPerPixelView<inner_view,func2_type>(dem_to_point_image( image.impl(), georef ), func2_type(georef.datum()) );
  }

} // end namespace vw

#endif//__DEMA_DEM_MANIPULATION_H__
