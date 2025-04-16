/*******************************************************************************************************************
								Link class to treat data from segmentation of images 2D (16bit only)
*******************************************************************************************************************/
#pragma once

typedef struct
{
	double x, y;
} val2D;

typedef struct
{
	size_t xMin, xMax, yMin, yMax;
} box2D;

typedef struct
{
	double maxPotentialWithNeighbor;
	int structureInContact;
	double minPotentialOnTheStructure;
	double maxPotentialOnTheStructure;
	int x;
	int y;
	bool kept;
	bool locked;
	int structureToBecaume;
	int newIndex;
} wsFusion2D;


typedef struct
{
	double val;
	size_t x, y;
} maxStruct2D;

typedef struct
{
	double xx, yy, xy;
} varianceStruct2D;

typedef struct
{
	double max, mean, angleMax;
} gaborCoefs;

typedef struct
{
	size_t peri;
	double surf;
	maxStruct2D min;
	maxStruct2D max;
	double stdev;
	double morpho;
	val2D moment;
	val2D momentSurf;
	double mean;
	double intensity;
	varianceStruct2D var;
	val2D sigma;
	double angle;
	size_t index;
	size_t time;
	int maxLength;
	size_t zOrderComputed;
	box2D boundingBox;
	size_t nbSubstructures;
} InfoDataRich2D;


typedef struct
{
	size_t peri;
	size_t surf;
	maxStruct2D max;
	double morpho;
	val2D moment;
	val2D momentSurf;
	double mean;
	double intensity;
	varianceStruct2D var;
	val2D sigma;
	double angle;
	size_t index;
	size_t time;
	int maxLength;
} InfoData2D;


class InfoStruct2D
{
public:
	explicit InfoStruct2D(const size_t index)
	{
		data_.index = index;
		data_.zOrderComputed = 0;
		next_ = nullptr;
	}

	~InfoStruct2D() = default;
	InfoStruct2D* next_;

	InfoDataRich2D data_;
};
