#ifndef TILING_H
#define TILING_H

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>


using namespace std;

struct TileConfig
{
	int tileWidth;
	int tileHeight;
	int xOverlap;
	int yOverlap;
	int featureMethod;
	int showResultsFlag;
	int tileScaleX;
	int tileScaleY;
};

class Tiling
{
public:
	//Constructor
	Tiling();
	
	//Destructor
	~Tiling();

	//Members
	TileConfig tileParams;
	IplImage* image;
	vector<CvRect> tiles;
	vector<IplImage*> tileImages;

	//Methods
	void printUsage();
	void printDefaultWarning(string configFilename);
	void readInputFiles(string inputFilename, vector<string>& inputFiles);
	void readTilingConfigFile(string filename);
	void saveTileBB();
	void process(IplImage* image);
	void setUpReferenceTiles(IplImage* image);
	void setUpMatchingTiles(Tiling ref, IplImage* im1);
	void clear();
};

#endif



