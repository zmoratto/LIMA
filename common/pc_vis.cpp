// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/visualization/pcl_visualizer.h>

int CountEntries(string pointCloudFilename)
{
    string line;
    ifstream myfile (pointCloudFilename.c_str());
    
    int numPoints = 0;
    if (myfile.is_open())
	{
    	while ( myfile.good() )
		{
			getline (myfile,line);
			float x = 0, y = 0, z = 0, r = 0, g = 0, b = 0;
			sscanf(line.c_str(), "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b);
			if ((x!=0.0)||(y!=0.0)||(z!=0.0))
			{
				numPoints++;
			}
      	}
    	myfile.close();
    }
    
    else cout << "Unable to open file." << endl;
    cout<<"numPoints="<<numPoints<<endl;
    return numPoints;
}

int main (int argc, char** argv)
{
	cout <<"num arguments = "<<argc<<endl;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
	string pointCloudFilename1, pointCloudFilename2;
	int numPoints;
	int subSampleStep = 1;
	int counter = 0;

	if (argc == 2)
	{
		pointCloudFilename1 = string(argv[1]);
		cout<<"pointCloudFilename1="<<pointCloudFilename1<<endl;
		numPoints = CountEntries(pointCloudFilename1);
	}
	if (argc == 3)
	{
		
		pointCloudFilename1 = string(argv[1]);
		cout<<"pointCloudFilename1="<<pointCloudFilename1<<endl;
		numPoints = CountEntries(pointCloudFilename1);
		subSampleStep = atoi(argv[2]);
	}

	if (argc == 1)
	{
		cout << endl << "Usage: ./pc_vis <pointcloudfilename> <step>" << endl << endl;
	}

	if ((argc == 2) || (argc==3))
	{
		cloud->width = numPoints/subSampleStep;
		cloud->height = 1;
		//cloud->is_dense = false;
		cloud->points.resize (cloud->width * cloud->height);

		if (numPoints>0)
		{
			ifstream myFile1 (pointCloudFilename1.c_str());
			string line;
			if (myFile1.is_open())
			{
				while ( myFile1.good() )
				{
					getline (myFile1,line);
					float x = 0, y = 0, z = 0, r = 255, g = 0, b = 0;
					sscanf(line.c_str(), "%f %f %f %f %f %f", &x, &y, &z, &r, &g, &b);
					if ((x!=0.0)||(y!=0.0)||(z!=0.0))
					{
						cloud->points[counter].x = x;
						cloud->points[counter].y = y;
						cloud->points[counter].z = z;
						cloud->points[counter].r = r;
						cloud->points[counter].g = g;
						cloud->points[counter].b = b;
						counter++;
					}
					for(int i=0; i<subSampleStep-1; i++)
						getline(myFile1, line);
				}
				myFile1.close(); 
			}
		}
	}
	else //Default demo of pc_vis
	{
		// Fill in the cloud data
		cloud->width= 5;
		cloud->height = 1;
		cloud->is_dense = false;
		cloud->points.resize (cloud->width * cloud->height);	
		cloud->points.resize (cloud->width * cloud->height);

		for (size_t i = 0; i < cloud->points.size (); ++i)
		{
			cloud->points[i].x = 1024 * rand () / (RAND_MAX + 1.0f);
			cloud->points[i].y = 1024 * rand () / (RAND_MAX + 1.0f);
			cloud->points[i].z = 1024 * rand () / (RAND_MAX + 1.0f);
		}

		for (size_t i = 0; i < cloud->points.size (); ++i)
			std::cerr << "" << cloud->points[i].x << " " << cloud->points[i].y << " " << cloud->points[i].z << std::endl;
	}

	cout << "Points Visualized: " << counter << endl;

	//display points
	pcl::visualization::CloudViewer viewer ("Simple Cloud Viewer");
	viewer.showCloud(cloud);
	while (!viewer.wasStopped()){}

	return (0);
}

