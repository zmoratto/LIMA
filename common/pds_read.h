// $Id$
// 
// pds2pgm: convert an unsigned short PDS image into a pgm image.
// 
// This program also reads from the input PDS image's header the CAHV
// camera parameters, translates them into a pan/tilt pair of values,
// and writes those into a comment in the output pgm image's header,
// that Viz can understand.
//
// Author: Doron Tal
// Date Created: Dec. 11, 2003

#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <limits.h>
#include <math.h>

extern "C" { 
#include "oal.h"
#include "lablib3.h"
}

using namespace std;

// #define _DEBUG

// location of the error report file created by the odl label parser 
#define ERRS_LOC "parser_errors.txt"

// only one bits/pixel value is supported for input images, defined here:
#define INPUT_BITS_PER_PIXEL 16

// Define pixels "outliers" to dump as potential noise
#define HISTOGRAM_OUTLIER_PERCENTAGE 0.005

// Use this to shift minimum brightness up a bit
#define BRIGHTNESS_OFFSET 0.0

// markers for string scanning
#define CAM_PARAMS_MARKER "GEOMETRIC_CAMERA_MODEL"
#define C_MARKER "MODEL_COMPONENT_1"
#define A_MARKER "MODEL_COMPONENT_2"
#define H_MARKER "MODEL_COMPONENT_3"
#define V_MARKER "MODEL_COMPONENT_4"
#define O_MARKER "MODEL_COMPONENT_5"
#define R_MARKER "MODEL_COMPONENT_6"

#define SITE_GEOM_MARKER "SITE_DERIVED_GEOMETRY_PARMS"
#define SUN_AZ_MARKER "SOLAR_AZIMUTH"
#define SUN_EL_MARKER "SOLAR_ELEVATION"

#define START_TIME_MARKER "SPACECRAFT_CLOCK_START_COUNT"
#define LR_MARKER "FRAME_ID"
#define RMC_MARKER "ROVER_MOTION_COUNTER"

#define ROVER_COORDINATE_MARKER "COORDINATE SYSTEM STATE: ROVER"
#define OFFSET_MARKER "ORIGIN_OFFSET_VECTOR"
#define ROTATION_MARKER "ORIGIN_ROTATION_QUATERNION"

//Left and right constants
#define RIGHT_ID 2
#define LEFT_ID 1
#define NEITHER_ID 0

int scan_file_for_string(FILE *fp, const char *str);
void a2pan_tilt(const double a[3], float &pan, float& tilt);

int scan_file_for_sun_angles(const char* strInFilename, 
			     float *azimuth, float *elevation);
int getRadianceDataForImage(ODLTREE odltree, double &radianceScaleFactor,
			    double &radianceOffset);
long get_pds_label_value(ODLTREE odltree, const char* strLabelName);
int scan_file_for_cahv_params(const char* strInFilename, 
                              float c[3], float a[3], float h[3], float v[3],
			      float o[3], float r[3], bool &gotOR);


int PDSReadImage(string strInFilename, unsigned char** pDataBuffer, 
                 long *r_nCols, long *r_nRows, long *r_nBands, long *r_nBits);

void Convert16BitTo8BitImage(unsigned short int* inputBuffer, unsigned char *outputBuffer, int nCols, int nRows);
void WriteCAHVORModel(string strCalFilename, float c[3], float a[3], float h[3], float v[3], float o[3], float r[3]);
void WriteCAHVModel(string strCalFilename, float c[3], float a[3], float h[3], float v[3]);
/*
bool write_pgm_byte_image(unsigned char *pBuffer, const char* filename, 
                          const long nCols, const long nRows, 
                          const float pan, const float tilt);
*/
int scan_file_for_clock_start(const char* strInFilename, double &start_time);

int scan_file_for_left_or_right(const char* strInFilename, int &result);
int scan_file_for_rover_PMA(const char* strInFilename, int &PMA);


int scan_file_for_origin_offset(const char* strInFilename, double offset[3]);
int scan_file_for_origin_rotation(const char* strInFilename, double quaternion[4]);

int scan_file_for_image_size( const char* strInFilename, long &nCols, long &nRows );


