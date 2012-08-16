#include <iostream>
#include <cmath>

using namespace std;


double
dotprod3(double v1[3], double v2[3])
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

void
crossprod3(double v1[3], double v2[3], double result[3])
{
  // a x b = (a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1)
  result[0] = (v1[1] * v2[2] - v1[2] * v2[1]);
  result[1] = (v1[2] * v2[0] - v1[0] * v2[2]);
  result[2] = (v1[0] * v2[1] - v1[1] * v2[0]);
}

double
magnitude(double v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

void
scalar_mult(double s, double v1[3], double result[3])
{
  result[0] = s * v1[0];
  result[1] = s * v1[1];
  result[2] = s * v1[2];
}

void
add_vec3(double v1[3], double v2[3], double result[3])
{
  result[0] = v1[0] + v2[0];
  result[1] = v1[1] + v2[1];
  result[2] = v1[2] + v2[2];
}

void
sub_vec3(double v1[3], double v2[3], double result[3])
{
  result[0] = v1[0] - v2[0];
  result[1] = v1[1] - v2[1];
  result[2] = v1[2] - v2[2];
}

void
PrintPinholeModel(double C[3], double A[3], double H[3], double V[3],
		  int cols, int rows, float pixelSize)
{
  double x_axis[3] = {0.0, 0.0, 0.0}, y_axis[3] = {0.0, 0.0, 0.0};
  double AxH[3], AxV[3], optical_center[2];
  crossprod3(A, H, AxH);
  crossprod3(A, V, AxV);
  double horizontal_scale = magnitude(AxH);
  double vertical_scale = magnitude(AxV);
  optical_center[0] = dotprod3(A, H);
  optical_center[1] = dotprod3(A, V);
  double fov_Y = 2.0 * atan2(rows, 2.0*vertical_scale);
  double fov_X = 2.0 * atan2(cols, 2.0*horizontal_scale);

  cout << "Horizontal scale = " << horizontal_scale << endl;
  cout << "Vertical scale = " << vertical_scale << endl;

  cout << "Focal length = " << ((horizontal_scale + vertical_scale)*0.5 *
				pixelSize)
       << endl;

  cout << "Optical center pixel = [" << lround(optical_center[0]) << ", "
       << lround(optical_center[1]) << "]" << endl;
  
  // Image plane x-axis: x_axis = (H - optical_center[0]*A)/horizontal_scale;
  scalar_mult(-optical_center[0], A, x_axis);
  add_vec3(H, x_axis, x_axis);
  scalar_mult(1.0/horizontal_scale, x_axis, x_axis);
  // Image plane y-axis: y_axis = (H - optical_center[1]*A)/vertical_scale;
  scalar_mult(-optical_center[1], A, y_axis);
  add_vec3(V, y_axis, y_axis);
  scalar_mult(1.0/vertical_scale, y_axis, y_axis);

  cout << "Center of projection = {"
       << C[0] << ", "
       << C[1] << ", "
       << C[2] << "}" << endl;

  cout << "Image plane X-axis = {"
       << x_axis[0] << ", "
       << x_axis[1] << ", "
       << x_axis[2] << "}" << endl;

  cout << "Image plane Y-axis = {"
       << y_axis[0] << ", "
       << y_axis[1] << ", "
       << y_axis[2] << "}" << endl;

  cout << "Image plane Z(A)-axis = {"
       << A[0] << ", "
       << A[1] << ", "
       << A[2] << "}" << endl;

  cout << "FOV X (deg) = " << fov_X * 180.0 / M_PI << endl;
  cout << "FOV Y (deg) = " << fov_Y * 180.0 / M_PI << endl;
}

void
GetPinholeModel(double C[3], double A[3], double H[3], double V[3],
		int cols, int rows, float pixelSize, double optical_center[2], double *fov_X, double *fov_Y)
{
  double x_axis[3] = {0.0, 0.0, 0.0}, y_axis[3] = {0.0, 0.0, 0.0};
  double AxH[3], AxV[3];//, optical_center[2];
  crossprod3(A, H, AxH);
  crossprod3(A, V, AxV);
  double horizontal_scale = magnitude(AxH);
  double vertical_scale = magnitude(AxV);
  optical_center[0] = dotprod3(A, H);
  optical_center[1] = dotprod3(A, V);
  /*double*/ *fov_Y = 2.0 * atan2(rows, 2.0*vertical_scale);
  /*double*/ *fov_X = 2.0 * atan2(cols, 2.0*horizontal_scale);
  /*
  cout << "Horizontal scale = " << horizontal_scale << endl;
  cout << "Vertical scale = " << vertical_scale << endl;

  cout << "Focal length = " << ((horizontal_scale + vertical_scale)*0.5 *
				pixelSize)
       << endl;

  cout << "Optical center pixel = [" << lround(optical_center[0]) << ", "
       << lround(optical_center[1]) << "]" << endl;
  */
  // Image plane x-axis: x_axis = (H - optical_center[0]*A)/horizontal_scale;
  scalar_mult(-optical_center[0], A, x_axis);
  add_vec3(H, x_axis, x_axis);
  scalar_mult(1.0/horizontal_scale, x_axis, x_axis);
  // Image plane y-axis: y_axis = (H - optical_center[1]*A)/vertical_scale;
  scalar_mult(-optical_center[1], A, y_axis);
  add_vec3(V, y_axis, y_axis);
  scalar_mult(1.0/vertical_scale, y_axis, y_axis);
  /*
  cout << "Center of projection = {"
       << C[0] << ", "
       << C[1] << ", "
       << C[2] << "}" << endl;

  cout << "Image plane X-axis = {"
       << x_axis[0] << ", "
       << x_axis[1] << ", "
       << x_axis[2] << "}" << endl;

  cout << "Image plane Y-axis = {"
       << y_axis[0] << ", "
       << y_axis[1] << ", "
       << y_axis[2] << "}" << endl;

  cout << "Image plane Z(A)-axis = {"
       << A[0] << ", "
       << A[1] << ", "
       << A[2] << "}" << endl;

  cout << "FOV X (deg) = " << *fov_X * 180.0 / M_PI << endl;
  cout << "FOV Y (deg) = " << *fov_Y * 180.0 / M_PI << endl;
  */
}


