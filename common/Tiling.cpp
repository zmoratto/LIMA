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
		cout << "*********************** USAGE **********************" << endl;
		cout << "./Tiling_test <configFile> <inputFile> <outputDir>  " << endl;
		cout << "****************************************************" << endl;
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
		fin >> identifier >> featureMethod;
		fin >> identifier >> showResultsFlag;
		fin >> identifier >> tileScaleX;
		fin >> identifier >> tileScaleY;

		fin.close();

		tileParams.tileWidth = tileWidth;
		tileParams.tileHeight = tileHeight;
		tileParams.xOverlap = xOverlap;
		tileParams.yOverlap = yOverlap;
		tileParams.featureMethod = featureMethod;
		tileParams.showResultsFlag = showResultsFlag;
		tileParams.tileScaleX = tileScaleX;
		tileParams.tileScaleY = tileScaleY;
	}
	else
	{
		tileParams.tileWidth = 500;
		tileParams.tileHeight = 500;
		tileParams.xOverlap = 50;
		tileParams.yOverlap = 50;
		tileParams.featureMethod = 0;
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

	for (r = 0; r < size.height; r += (y-yOver))
	{
		for (c = 0; c < size.width; c += (x-xOver))
		{
			roi.x = c;
			roi.y = r;
			if(c+x>size.width)
				roi.width = size.width-c;
			else
				roi.width = x;
			if(r+y>size.height)
				roi.height = size.height-r;
			else
				roi.height = y;

			if(roi.width <= xOver)
				break;
			if(roi.height <= yOver)
				break;

			tiles.push_back(roi);

			if(x >= size.width)
				break;
		}

		if(y >= size.height)
			break;
	}
}

void Tiling::setUpMatchingTiles(Tiling ref, IplImage* im1)
{
	int xScale = ref.tileParams.tileScaleX;
	int yScale = ref.tileParams.tileScaleY;
	CvRect temp;

	for(int i=0;i<ref.tiles.size();i++)
	{
		if(ref.tiles[i].x - xScale < 0)
		{
			temp.width = ref.tiles[i].width + ref.tiles[i].x;
			temp.x = 0;
		}
		else
		{
			temp.width = ref.tiles[i].width + xScale;
			temp.x = ref.tiles[i].x - xScale;
		}

		if(ref.tiles[i].y - yScale <0)
		{
			temp.height = ref.tiles[i].height + ref.tiles[i].y;
			temp.y = 0;
		}
		else
		{
			temp.height = ref.tiles[i].height + yScale;
			temp.y = ref.tiles[i].y - yScale;
		}

		tiles.push_back(temp);

		if(tiles[i].x + tiles[i].width + xScale > im1->width)
		{
			tiles[i].width = im1->width - tiles[i].x;
		}
		else
		{
			tiles[i].width += xScale;
		}

		if(tiles[i].y + tiles[i].height + yScale > im1->height)
		{
			tiles[i].height = im1->height - tiles[i].y;
		}
		else
		{
			tiles[i].height += yScale;
		}
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
