#include <iostream>

#include <Camera.h>
#include <Pvl.h>
#include <CameraFactory.h>
#include <CameraDetectorMap.h>

using namespace std;

vector<float> map2cam(std::string map_projected_cube, std::vector<float> map_pixel)
{
 
  std::vector<float> cam_pixel;
  cam_pixel.resize(2);

  // Load up ISIS objects about the map projected cube
  Isis::Pvl label( map_projected_cube );
  Isis::Camera* camera = Isis::CameraFactory::Create( label );

  // Note that ISIS is different from C. They start their index at 1.

  // convert a map projected pixel location to the
  // original image coordinate system.
  camera->SetImage(map_pixel[0], map_pixel[1]);
  cam_pixel[0] = camera->DetectorMap()->ParentSample();
  cam_pixel[1] = camera->DetectorMap()->ParentLine();

  // delete remaining ISIS objects
  delete camera;

  return cam_pixel;
}
