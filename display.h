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

#endif /* DISPLAY_H */
