// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
#include <vw/Math/Vector.h>
#include <vw/Math/Matrix.h>

#include <stdio.h>
#include <asp/IsisIO.h>
#include <asp/IsisIO/IsisCameraModel.h>

#include <Camera.h>
#include <Pvl.h>
#include <CameraFactory.h>

using namespace vw;
using namespace vw::math;
using namespace vw::cartography;

using namespace std;
#include <math.h>
#include "coregister.h"
#include "tracks.h"
#include "match.h"


//displays the original tracks(red) and transformed tracks(cyan) over the image
void ShowFinalTrackPtsOnImage(vector<vector<LOLAShot> >trackPts, Vector<float, 6> d, 
                              vector<int> trackIndices, string cubFilename, string outFilename)
{

  boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
  double nodataVal = rsrc->nodata_read();
  cout<<"nodaval:"<<nodataVal<<endl;
  DiskImageView<PixelGray<float> > cub( rsrc );
  camera::IsisCameraModel model(cubFilename);
 

  vector<pointCloud> ptHere;
  int minX = cub.cols()-1;
  int maxX = 0;
  int minY = cub.rows()-1;
  int maxY = 0;
  int point_size = 7;

  //for (int i = 0; i < trackPts.size(); i++){//for each track
 for (unsigned int i = 0; i < trackIndices.size(); i++){//for each selected track 
    for(unsigned int j = 0; j < trackPts[i].size(); j++){ //for each shot in a track
      for(unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){ //for each pt in a shot 
 
	  pointCloud pt = trackPts[i][j].LOLAPt[k]; 
	  float lon = pt.x();
	  float lat = pt.y();
	  float rad = pt.z();
       
          Vector3 lon_lat_rad (lon,lat,rad*1000);
          Vector3 xyz = lon_lat_radius_to_xyz(lon_lat_rad);
          Vector2 DRG_pix = model.point_to_pixel(xyz);
          float x = DRG_pix[0];
          float y = DRG_pix[1];

          //compute the transformed pts
          int xd = (int)floor(d[0]*x + d[1]*y + d[2]);
          int yd = (int)floor(d[3]*x + d[4]*y + d[5]);

          if (x < minX){
		minX = x;
          }
	  if (y < minY){
		minY = y;
	  }
	  if (x > maxX){
		maxX = x;
	  }
	  if (y > maxY){
		maxY = y;
	  }
	  
	  if (xd < minX){
		minX = xd;
	  }
	  if (yd < minY){
		   minY = yd;
	  }
	  if (xd > maxX){
		maxX = xd;
           }
	  if (yd > maxY){
		maxY = yd;
	  }

        
	}
      }
    }
   
   //extend the bounding box
   maxX = maxX + 4*point_size;
   minX = minX - 4*point_size;
   maxY = maxY + 4*point_size;
   minY = minY - 4*point_size;

   //make sure the bounding box is within the image boundaries
   if (minX < 0) minX = 0;   
   if (minY < 0) minY = 0; 
   if (maxX > cub.cols()-1) maxX = cub.cols()-1;   
   if (maxY > cub.rows()-1) maxY = cub.rows()-1; 
  
   printf("minX = %d, minY = %d, maxX = %d,  maxY = %d\n", minX, minY, maxX, maxY);

   ImageView<PixelGray<uint8> > DRG = apply_mask(normalize(create_mask(cub,nodataVal))*255,0);
   ImageView<PixelRGB<uint8> > DRG_crop = crop(DRG, int32(minX), int32(minY), maxX-minX+1, maxY-minY+1);
   
   //#if 0
   for (unsigned int i = 0; i < trackPts.size(); i++){//for each track
    for(unsigned int j = 0; j < trackPts[i].size(); j++){ //for each shot in a track
      for(unsigned int k = 0; k < trackPts[i][j].LOLAPt.size(); k++){ //for each pt in a shot 
       
	if ((trackPts[i][j].valid == 1) && (trackPts[i][j].featurePtLOLA == 1)){ 
            int xl, yt, w, h;

		pointCloud pt = trackPts[i][j].LOLAPt[k]; 
		float lon = pt.x();
		float lat = pt.y();
		float rad = pt.z();
	  
            Vector3 lon_lat_rad (lon,lat,rad*1000);
            Vector3 xyz = lon_lat_radius_to_xyz(lon_lat_rad);
            Vector2 DRG_pix = model.point_to_pixel(xyz);
            float x = DRG_pix[0];
            float y = DRG_pix[1];
                  
            //make sure we are not running outside the image boundaries.
            xl = int32(x) - point_size-minX;
            yt = int32(y) - point_size-minY;
            w = point_size;
            h = point_size;
          
            if ((xl>0) && (yt>0) && (xl<DRG_crop.cols()-point_size-1) && (yt<DRG_crop.rows()-point_size-1)){
                fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(255, 0, 0));
		}
        
		//compute the transformed pts
		int xd = (int)floor(d[0]*x + d[1]*y + d[2]);
		int yd = (int)floor(d[3]*x + d[4]*y + d[5]);

            //make sure the matched point is inside the image boundaries.
		xl = int32(xd) - point_size-minX;
            yt = int32(yd) - point_size-minY;
            w = point_size;
            h = point_size;

            if ((xl>0) && (yt>0) && (xl<DRG_crop.cols()-point_size-1) && (yt<DRG_crop.rows()-point_size-1)){
               fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(0, 0, 255));

              /*  
              if (trackPts[i][j].featurePtRefl == 1.0){
		fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(0, 255, 0));
              }
		  else{
		fill(crop(DRG_crop, xl, yt, w, h), PixelRGB<uint8>(0, 0, 255));
		  }
		  */
		}
      
	}//valid track && LOLA feature point
      }//each point in shot
   
    }//each shot
   }
   //#endif
  //write output image
  GeoReference moonref( Datum("D_MOON"), identity_matrix<3>() );
   
  write_georeferenced_image(outFilename,
                            DRG_crop,
                            moonref, TerminalProgressCallback("Core","Processing:"));
   
}

void SaveBigGCPImages(vector<gcp> gcps,  string cubFile, string filename)
{
	if (gcps.size() <= 0)
		return;
	cout << "Saving file " << filename << ".\n";
	boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFile));
	double nodataVal = rsrc->nodata_read();
	DiskImageView<PixelGray<float> > cub( rsrc );
	int width = cub.cols();
	int height = cub.rows();
	
	int point_size = 7;
	int w = point_size;
	int h = point_size;
	int xl, yt;
	
	ImageView<PixelRGB<uint8> > assembledImg(width, height);
	assembledImg = apply_mask(normalize(create_mask(cub,nodataVal))*255,0);
	
	for (unsigned int j = 0; j < gcps.size(); j++)
	{
		unsigned int i;
		for (i = 0; i < gcps[j].filename.size(); i++)
				if (gcps[j].filename[i] == cubFile)
						break;
		if (i >= gcps[j].filename.size())
			continue;
		//draw the interest point before alignment
		xl = int32(gcps[j].x_before[i]) - point_size / 2;
		yt = int32(gcps[j].y_before[i]) - point_size / 2;
		if ((xl>0) && (yt>0) && (xl < assembledImg.cols()-point_size-1) && (yt < assembledImg.rows()-point_size-1))
			fill(crop(assembledImg, xl, yt, w, h), PixelRGB<uint8>(255, 0, 0));
		
		//draw the interest point after alignment
		xl = int32(gcps[j].x[i]) - point_size;
		yt = int32(gcps[j].y[i]) - point_size;
		if ((xl>0) && (yt>0) && (xl < assembledImg.cols()-point_size-1) && (yt < assembledImg.rows()-point_size-1))
			fill(crop(assembledImg, xl, yt, w, h), PixelRGB<uint8>(0, 0, 255));
	}
	
	//save the assembled image to file
	write_image(filename, assembledImg);
}

//displays info in GCP file
void SaveGCPImages(struct gcp this_gcp, string assembledImgFilename)
{
 

  //create the assembled image;
  int numHorBlocks = 4;
  int numVerBlocks = 4;
  int index = 0;
  int colIndex = 0;
  int rowIndex = 0;
  int blockWidth = 200;
  int blockHeight = 200;
  int point_size = 7;

  numVerBlocks = (this_gcp.filename.size()/numHorBlocks);
  if ((unsigned int)(numVerBlocks*numHorBlocks) < this_gcp.filename.size()){
    numVerBlocks = numVerBlocks + 1;
  }

  ImageView<PixelRGB<uint8> > assembledImg(numHorBlocks*blockWidth, numVerBlocks*blockHeight);

  for (unsigned int i = 0; i < this_gcp.filename.size(); i++){
    
    string cubFilename = this_gcp.filename[i];
    boost::shared_ptr<DiskImageResource> rsrc( new DiskImageResourceIsis(cubFilename) );
    double nodataVal = rsrc->nodata_read();
    cout<<"nodaval:"<<nodataVal<<endl;
    DiskImageView<PixelGray<float> > cub( rsrc );
    int width = cub.cols();
    int height = cub.rows();
  
 
    cout<<"before_x="<<this_gcp.x_before[i]<<" before_y="<<this_gcp.y_before[i]<<endl;
    cout<<"after_x="<<this_gcp.x[i]<<" after_y="<<this_gcp.y[i]<<endl;

    float minX, minY, maxX, maxY;
    /*
    //centered around the aligned features
    minX = xPosArray[i]-blockWidth/2;
    maxX = xPosArray[i]+blockWidth/2;
    minY = yPosArray[i]-blockHeight/2;
    maxY = yPosArray[i]+blockHeight/2;
    */
   
    //centered around the original features
    minX = this_gcp.x_before[i]-blockWidth/2;
    maxX = this_gcp.x_before[i]+blockWidth/2;
    minY = this_gcp.y_before[i]-blockHeight/2;
    maxY = this_gcp.y_before[i]+blockHeight/2;

    int adjustedBlockWidth = blockWidth;
    int adjustedBlockHeight = blockHeight; 
    int w = point_size;
    int h = point_size;
    int xl, yt;

    if (maxY > height-1){adjustedBlockHeight = blockHeight-(maxY-height+1);}   
    if (maxX > width-1){adjustedBlockWidth = blockWidth - (maxX-width+1);}
    if (minY < 0){minY = 0;}
    if (minX < 0){minX = 0;}
    
    if ((adjustedBlockHeight > 0) && (adjustedBlockWidth > 0)){
      
      crop( assembledImg, colIndex*blockWidth, rowIndex*blockHeight, adjustedBlockWidth, adjustedBlockHeight ) = 
      crop(apply_mask(normalize(create_mask(cub,nodataVal))*255,0), int32(minX), int32(minY), adjustedBlockWidth, adjustedBlockHeight); 
    
      //draw the interest point before alignment
      xl = int32(this_gcp.x_before[i]) - minX - point_size + colIndex*blockWidth;
      yt = int32(this_gcp.y_before[i]) - minY - point_size + rowIndex*blockHeight;    
      if ((xl>0) && (yt>0) && (xl < assembledImg.cols()-point_size-1) && (yt < assembledImg.rows()-point_size-1)){
        fill(crop(assembledImg, xl, yt, w, h), PixelRGB<uint8>(255, 0, 0));
      }

      //draw the interest point after alignment
      xl = int32(this_gcp.x[i]) - minX - point_size + colIndex*blockWidth;
      yt = int32(this_gcp.y[i]) - minY - point_size + rowIndex*blockHeight;
      if ((xl>0) && (yt>0) && (xl < assembledImg.cols()-point_size-1) && (yt < assembledImg.rows()-point_size-1)){
        fill(crop(assembledImg, xl, yt, w, h), PixelRGB<uint8>(0, 0, 255));
      }

    }

    index++;
    rowIndex = index/numHorBlocks;
    colIndex = index - rowIndex*numHorBlocks;
     
   
  }

  //save the assembled image to file
  write_image(assembledImgFilename, assembledImg);

}



//displays the LOLA tracks in an image format on a black background
//quick space efficient way to visualize Lidar data without storing the real image but a black background. 
void MakeGrid(vector<vector<LOLAShot> >trackPts, int numVerPts, int numHorPts, string DEMFilename, vector<int> trackIndices)
{
  int l, m/*, n*/;  
  ImageView<PixelGray<float> > DEMImage(numHorPts, numVerPts);
  GeoReference DEMgeo;

  
  if (trackIndices.size() == 1){
      //string outDEMFilename;
      //char* outDEMFilename_char = new char[500];
      //sprintf (outDEMFilename_char, "../results/dem_%d.tiff", k);
      //outDEMFilename = std::string(outDEMFilename_char);
  }

  Vector4 coords = FindMinMaxLat(trackPts);

  printf("MakeGrid: minLat=%f, maxLat=%f, minLon=%f maxLon=%f\n", coords(0), coords(1), coords(2), coords(3));

  float minLat = coords(0); //this causes seg-fault
  float maxLat = coords(1); 
  float minLon = coords(2); 
  float maxLon = coords(3);

  float lonDelta = (maxLon-minLon)/numHorPts;
  float latDelta = (maxLat-minLat)/numVerPts;

  for (l = 0; l < numHorPts; l++){
    for (m = 0; m < numVerPts; m++){
      DEMImage(l, m) = 0.0;
    }
  }
  //fill the DEM

  printf("MakeGrid: numTracks = %d\n", (int)(trackIndices.size()));
  for (unsigned int k = 0; k < trackIndices.size();k++){
    int trackIndex = trackIndices[k];
    printf("MakeGrid: trackIndex = %d\n", trackIndex);
    for (unsigned n = 0; n < trackPts[trackIndex].size(); n++){ 
      for (unsigned int s = 0; s < trackPts[trackIndex][n].LOLAPt.size(); s++){

        float lon_index = (trackPts[trackIndex][n].LOLAPt[s].x() - minLon)/lonDelta;
        float lat_index = (trackPts[trackIndex][n].LOLAPt[s].y() - minLat)/latDelta;
        l = (int)floor(lon_index);
        m = (int)floor(lat_index);

        if ((m < numVerPts) && (l< numHorPts)){ 
          DEMImage(l, m) = trackPts[trackIndex][n].LOLAPt[s].z(); 
        }
        else{
          printf("Make Grid: Error\n");
          printf("Make Grid: l = %d, m = %d, numHorPts = %d, numVerPts = %d\n", l, m, numHorPts, numVerPts);
        }
      }
    }
  }
  write_georeferenced_image(DEMFilename, DEMImage,
                            DEMgeo, TerminalProgressCallback("{Core}","Processing:"));

}

DiskImageResource* get_image(string file, bool is_cub)
{
	if (is_cub)
		return new DiskImageResourceIsis(file);
	else
		return new DiskImageResourceGDAL(file);
}


/*void SaveAdjustedReflectanceImages(vector<gcp> gcps, vector<vector<LOLAShot> > tracks,  string cubFile, string filename, int image_id, bool is_cub)
{
	cout << "Saving adjusted reflectance image " << filename << ".\n";
	DiskImageResource* image =  get_image(cubFile, is_cub);
	boost::shared_ptr<DiskImageResource> rsrc((DiskImageResource*)(image));
	DiskImageView<PixelGray<float> > cub(rsrc);
	
	int width = cub.cols();
	int height = cub.rows();
	ImageView<PixelRGB<uint8> > assembledImg(width, height);
	
	if (is_cub)
	{
		double nodataVal = rsrc->nodata_read();
		assembledImg = apply_mask(normalize(create_mask(cub,nodataVal))*255,0);
	}
	// calibrated image
	else
	{
		assembledImg = apply_mask(normalize(cub) * 255, 0);
	}

	int xl, yt;
	
	unsigned int gcp_index = 0;
	for (unsigned int i = 0; i < tracks.size(); i++)
	{
		int last_x = -1;
		int next_x = -1;
		int last_y = -1;
		int next_y = -1;
		int last_y_before = -1;
		int next_y_before = -1;
		while (gcp_index < gcps.size() && (gcps[gcp_index].x[image_id] < 0)) gcp_index++;
		if (gcp_index < gcps.size() && (unsigned int)gcps[gcp_index].trackIndex == i)
		{
			next_x = gcps[gcp_index].x[image_id];
			next_y = gcps[gcp_index].y[image_id];
			next_y_before = gcps[gcp_index].y_before[image_id];
		}
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			LOLAShot* s = &tracks[i][j];
			if (!s->valid)
				continue;
			xl = s->image_x;
			yt = s->image_y;

			bool key_point = false;
			if (gcp_index < gcps.size() && (unsigned int)gcps[gcp_index].trackIndex == i && (unsigned int)gcps[gcp_index].shotIndex == j)
			{
				do gcp_index++; while(gcp_index < gcps.size() && gcps[gcp_index].x[image_id] < 0);
				if (gcp_index < gcps.size() && (unsigned int)gcps[gcp_index].trackIndex == i)
				{
					last_x = next_x;
					next_x = gcps[gcp_index].x[image_id];
					last_y = next_y;
					next_y = gcps[gcp_index].y[image_id];
					last_y_before = next_y_before;
					next_y_before = gcps[gcp_index].y_before[image_id];
					key_point = true;
				}
				else
				{
					last_x = -1;
					next_x = -1;
					last_y = -1;
					next_y = -1;
					next_y_before = -1;
					last_y_before = -1;
				}
			}
			if (key_point && xl>0 && yt>0 && (xl < assembledImg.cols()-5) && (yt < assembledImg.rows()-2))
			{
				fill(crop(assembledImg, xl-3, yt, 11, 2), PixelRGB<uint8>(255, 0, 0));
				uint8 col = 255 * (s->luminence);
				fill(crop(assembledImg, xl+1, yt, 3, 2), PixelRGB<uint8>(col, col, col));
			}

			if (last_x != -1)
			{
				xl = int32(last_x + ((float)(yt - last_y_before)) / (next_y_before - last_y_before) * (next_x - last_x)) - 2;
				yt = int32(last_y + ((float)(yt - last_y_before)) / (next_y_before - last_y_before) * (next_y - last_y)) - 1;
			}
    
			if (xl>0 && yt>0 && (xl < assembledImg.cols()-5) && (yt < assembledImg.rows()-2))
			{
				PixelRGB<uint8> col(0, 0, 255);
				if (key_point)
				{
					col = PixelRGB<uint8>(0, 255, 0);
					fill(crop(assembledImg, xl-5, yt, 15, 2), col);
					fill(crop(assembledImg, xl+1, yt, 3, 2), PixelRGB<uint8>(uint8(255 * s->reflectance), uint8(255 * s->reflectance), uint8(255 * s->reflectance)));
				}
				else
				{
					fill(crop(assembledImg, xl, yt, 5, 2), col);
					fill(crop(assembledImg, xl+1, yt, 3, 2), PixelRGB<uint8>(uint8(255 * s->reflectance), uint8(255 * s->reflectance), uint8(255 * s->reflectance)));
				}
			}
			
		}
	}
	
	//save the assembled image to file
	write_image(filename, assembledImg);
}*/

void SaveReflectanceImages(vector<vector<AlignedLOLAShot> >& tracks,  ImageView<PixelGray<float> > cub, string filename)
{
	int width = cub.cols();
	int height = cub.rows();
	ImageView<PixelRGB<uint8> > assembledImg(width, height);
	
	assembledImg = apply_mask(normalize(cub) * 255, 0);

	int xl, yt;
	
	for (unsigned int i = 0; i < tracks.size(); i++)
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			AlignedLOLAShot* s = &tracks[i][j];
			if (s->synth_image < 0 || s->image < 0 || !s->valid)
				continue;
			xl = s->image_x;
			yt = s->image_y;
    
			if (xl>0 && yt>0 && (xl < assembledImg.cols()-5) && (yt < assembledImg.rows()-2))
			{
				fill(crop(assembledImg, xl, yt, 5, 2), PixelRGB<uint8>(255, 0, 0));
				uint8 col = (uint8)(255 * s->synth_image);
				//printf("%g %g %g\n", s->reflectance, t, cub(xl, yt)[0]);
				fill(crop(assembledImg, xl+1, yt, 3, 2), PixelRGB<uint8>(col, col, col));
			}
		}
	
	//save the assembled image to file
	write_image(filename, assembledImg);
}

void Save3DImage(vector<vector<LOLAShot> >& tracks, string filename)
{
	char buf[100];

	cout << "Saving 3D image file " << filename << ".\n";
	
	FILE* fv = fopen("vertices.tmp", "w");
	FILE* ff = fopen("faces.tmp", "w");
	unsigned int vertex = 1;

	for (unsigned int i = 0; i < tracks.size(); i++)
		for (unsigned int j = 0; j < tracks[i].size(); j++)
		{
			LOLAShot* s = &tracks[i][j];
			unsigned int num = s->LOLAPt.size();
			if (num < 3)
				continue;
			for (unsigned int k = 0; k < num; k++)
			{
        			Vector3 xyz = lon_lat_radius_to_xyz(s->LOLAPt[k]);
				fprintf(fv, "v %g %g %g\n", xyz.x(), xyz.y(), xyz.z());
			}
			for (unsigned int k = vertex; k < vertex + num - 2; k++)
			{
				fprintf(ff, "f %d %d %d\n", k, k+1, k+2);
				fprintf(ff, "f %d %d %d\n", k+2, k+1, k);
			}
			fprintf(ff, "f %d %d %d\n", vertex + num-2, vertex + num-1, vertex);
			fprintf(ff, "f %d %d %d\n", vertex, vertex + num-1, vertex + num - 2);
			fprintf(ff, "f %d %d %d\n", vertex + 1, vertex, vertex + num - 1);
			vertex += num;
		}
	snprintf(buf, 100, "cat vertices.tmp faces.tmp > %s", filename.c_str());
	std::stringstream out2;
	out2 << "cat vertices.tmp faces.tmp > " << filename;
	fclose(fv);
	fclose(ff);
	int ret = system(out2.str().c_str());
	ret = ret || system("rm vertices.tmp");
	ret = ret || system("rm faces.tmp");
	if (ret)
		fprintf(stderr, "Error saving vertex file.\n");
}

void overlay_image(char* image1, char* image2, Matrix3x3 H, char* outFile)
{
	boost::shared_ptr<DiskImageResource> rsrc1(new DiskImageResourceIsis(image1));
	boost::shared_ptr<DiskImageResource> rsrc2(new DiskImageResourceIsis(image2));
	DiskImageView<PixelGray<float> > temp1(rsrc1), temp2(rsrc2);

	ImageView<PixelGray<float> > base = normalize(temp1);
	ImageView<PixelGray<float> > t = normalize(temp2);

	for (int i = 0; i < base.rows(); i++)
		for (int j = 0; j < base.cols(); j++)
		{
			Vector3 coord = H * Vector3(i, j, 1);
			coord = coord / coord(2);
			int x = (int)coord(0);
			int y = (int)coord(1);
			if (x < 0 || y < 0 || x >= t.rows() || y >= t.cols())
				continue;
			base(i, j) = 0.5 * base(i, j) + 0.5 * t(x, y);
		}
	
	ImageView<PixelRGB<uint8> > assembledImg(temp1.rows(), temp1.cols());
	assembledImg = apply_mask(base * 255, 0);
	write_image(outFile, assembledImg);
}

