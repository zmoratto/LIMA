#include "Tiling.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

Tiling::Tiling()
{
}

Tiling::~Tiling()
{
}

void Tiling::printUsage()
{
		cout << endl;
		cout << "****************** USAGE *****************" << endl;
		cout << "./Tiling_test <configFile> <inputFile>              " << endl;
		cout << "******************************************" << endl;
		cout << endl;
}

void Tiling::printDefaultWarning(string configFilename)
{
		cout << endl;
		cout << "*************************************************************" << endl;
		cout << "WARNING: Unable to open file " << configFilename << "..." << endl;
		cout << "Using default parameters..." << endl;
		cout << "*************************************************************" << endl;
		cout << endl;
}

void Tiling::readInputFiles(string inputFilename, vector<string>& inputFiles)
{
	ifstream fin;
	string line;

	fin.open(inputFilename.c_str());

	while(fin.good())
	{
		getline(fin,line);
		inputFiles.push_back(line);		
	}
	fin.close();
}

void Tiling::readTilingConfigFile(string filename)
{
	string line;
	ifstream fin;
	string identifier;
	int tileWidth;
	int tileHeight;
	int xOverlap;
	int yOverlap;
	int featureMethod;
	int showResultsFlag;
	int tileScaleX, tileScaleY;
	fin.open(filename.c_str());

	if(fin.is_open())
	{
		fin >> identifier >> tileWidth;
		fin >> identifier >> tileHeight;
		fin >> identifier >> xOverlap;
		fin >> identifier >> yOverlap;


		fin.close();

		tileParams.tileWidth = tileWidth;
		tileParams.tileHeight = tileHeight;
		tileParams.xOverlap = xOverlap;
		tileParams.yOverlap = yOverlap;
	}
	else
	{
		tileParams.tileWidth = 500;
		tileParams.tileHeight = 500;
		tileParams.xOverlap = 50;
		tileParams.yOverlap = 50;
		printDefaultWarning(filename);
	}
}

void Tiling::saveTileBB()
{
	string filename;
	ofstream fout;
	filename = string("Tile_BB.txt");
	fout.open(filename.c_str());
	for(int i=0; i<tiles.size(); i++)
	{
		fout << i << "\t" << tiles[i].x << " " << tiles[i].y << " " << tiles[i].width << " " << tiles[i].height << endl;
	}
	fout.close();
}


void Tiling::process(IplImage* im1)
{
	tileImages.clear();

	for(int i=0;i<tiles.size();i++)
	{
		IplImage* im2 = cvCloneImage(im1);
	    cvSetImageROI(im2, tiles[i]);
		IplImage* copy = cvCloneImage(im2);
		tileImages.push_back(copy);
		tileImages[i]->width = tiles[i].width;
		tileImages[i]->height = tiles[i].height;
		cvReleaseImage(&im2);
	}
}

void Tiling::setUpReferenceTiles(IplImage* image)
{
	CvRect roi;
	CvSize size;
	int r,c,x,y,oldC,oldR;
	int xOver, yOver;

	x = tileParams.tileWidth;
	y = tileParams.tileHeight;
	xOver = tileParams.xOverlap;
	yOver = tileParams.yOverlap;

	size.height = image->height;
	size.width = image->width;

	tiles.clear();

   //Check if Tile too large
   if( x > size.width || x < 0)
      x = size.width;
   if( y > size.height || y < 0 )
      y = size.height;

	for (r = 0; r < size.height; r += (y-yOver))
    {
        for (c = 0; c < size.width; c += (x-xOver))
        {
         roi.width = x;
         roi.height = y;
         if(c+x>size.width)
             roi.x = size.width-x;
         else
             roi.x = c;
         if(r+y>size.height)
             roi.y = size.height-y;
         else
             roi.y = r;

         tiles.push_back(roi);

         if( c+x >= size.width )
            break;

        }

    if( r+y >= size.height )
         break;

    }

}

void Tiling::clear()
{
	for(int i=0; i<tileImages.size(); i++)
	{
		cvReleaseImage(&tileImages[i]);
	}

	tileImages.clear();
}
