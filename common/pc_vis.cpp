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
			//cout << line << endl;
			float x, y, z, g;
			sscanf(line.c_str(), "%f %f %f", &x, &y, &z);
			//sscanf(line.c_str(), "%f %f %f %f", &x, &y, &z &g);
			if ((x!=0.0)||(y!=0.0)||(z!=0.0))
			{
				//cout<<"x="<<x<<", y="<<y<<", z="<<z<<endl;
				numPoints++;
			}
      	}
    	myfile.close();
    }
    
    else cout << "Unable to open file";
    cout<<"numPoints="<<numPoints<<endl;
    return numPoints;

}

int main (int argc, char** argv)
{
	cout <<"num arguments = "<<argc<<endl;
	//pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGB>);
	string pointCloudFilename1, pointCloudFilename2;
	int numPoints, numPoints1, numPoints2;
	numPoints1 = 0;
	numPoints2 = 0;

	if (argc == 2)
	{
		pointCloudFilename1 = string(argv[1]);
		cout<<"pointCloudFilename1="<<pointCloudFilename1<<endl;
		numPoints1 = CountEntries(pointCloudFilename1);
	}
	if (argc == 3)
	{
		pointCloudFilename1 = string(argv[1]);
		cout<<"pointCloudFilename1="<<pointCloudFilename1<<endl;
		numPoints1 = CountEntries(pointCloudFilename1);

		pointCloudFilename2 = string(argv[2]);
		cout<<"pointCloudFilename2="<<pointCloudFilename2<<endl;
		numPoints2 = CountEntries(pointCloudFilename2);
	}

	if (argc == 1)
	{
		cout << endl << "This program expects at least one point cloud file name as a command line argument." << endl << endl;
		return -1;
	}

	numPoints = numPoints1 + numPoints2;
	 
	if ((argc == 2) || (argc==3))
	{
		cloud->width = numPoints;
		cloud->height = 1;
		//cloud->is_dense = false;
		cloud->points.resize (cloud->width * cloud->height);

		int counter = 0;

		if (numPoints1>0)
		{
			ifstream myFile1 (pointCloudFilename1.c_str());
			string line;
			if (myFile1.is_open())
			{
				while ( myFile1.good() )
				{
					getline (myFile1,line);
					float x, y, z, r, g, b;
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
				}
				myFile1.close(); 
			}
		}
	}
	else
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

		//pcl::io::savePCDFileASCII ("test_pcd.pcd", *cloud);
		//std::cerr << "Saved " << cloud->points.size () << " data points to test_pcd.pcd." << std::endl;

		for (size_t i = 0; i < cloud->points.size (); ++i)
			std::cerr << "" << cloud->points[i].x << " " << cloud->points[i].y << " " << cloud->points[i].z << std::endl;
	}

	//display points
	pcl::visualization::CloudViewer viewer ("Simple Cloud Viewer");
	viewer.showCloud(cloud);
	while (!viewer.wasStopped()){}

	return (0);
}

