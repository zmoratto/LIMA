// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef DISPLAY_H
#define DISPLAY_H

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;

//these functions will be merged into tracks class-START
void ShowFinalTrackPtsOnImage(vector<vector<LOLAShot> >trackPts, Vector<float, 6> d, 
                              vector<int> trackIndices, string DRGFilename, string outFilename);
void MakeGrid(vector<vector<LOLAShot> >trackPts, int numVerPts, int numHorPts, string DEMFilename, vector<int> trackIndices);
//these functions will be merged into tracks class-END

void SaveGCPImages(string GCPFilename, string cubDirname, string assembledImgFilename);
void SaveGCPImages(struct gcp this_gcp,  string assembledImgFilename);
void SaveBigGCPImages(vector<gcp> gcps, string cubFile, string assembledImgFilename);

/**
 * Draw aligned tracks on top of image, and save to file imgFilename.
 **/
void SaveReflectanceImages(vector<vector<AlignedLOLAShot> >& trackPts,  ImageView<PixelGray<float> > cubFile, string imgFilename);
/**
 * Save a mosaic of two overlapping images given a transformation from coordinates
 * in image1 to image2.
 **/
void overlay_image(char* image1, char* image2, Matrix3x3 H, char* outFile);

#endif /* DISPLAY_H */
