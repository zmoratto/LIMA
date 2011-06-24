#include <iostream>

// Bring in ISIS headers
#include <Camera.h>
#include <Pvl.h>
#include <CameraFactory.h>
#include <CameraDetectorMap.h>
using namespace std;

vector<float> map2cam(std::string map_projected_cube, std::vector<float> map_pixel);

