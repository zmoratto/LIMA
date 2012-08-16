#include <iostream>
#include <cmath>

void PrintPinholeModel(double C[3], double A[3], double H[3], double V[3],
		  int cols, int rows, float pixelSize);
void
GetPinholeModel(double C[3], double A[3], double H[3], double V[3],
		int cols, int rows, float pixelSize, double optical_center[2], double *fov_X, double *fov_Y);
