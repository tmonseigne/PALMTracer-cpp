#ifndef __link_2D__
#define __link_2D__

#include <stdio.h>
#include <math.h>
#include "allDataLib.h"

extern int __nbInfoStrucAlloc;
extern int __nbInfoStrucDealloc;


#define MAX_NB_ANGLES_GABOR	8
#define MAX_NB_FREQ_GABOR	13

typedef struct
{
	double x;
	double y;
} val2D;


typedef struct
{
	long x;
	long y;
} coor2D;

typedef struct
{
	long xMin;
	long xMax;
	long yMin;
	long yMax;
} box2D;

typedef struct
{
	float maxPotentialWithNeighbor;
	long structureInContact;
	float minPotentialOnTheStructure;
	float maxPotentialOnTheStructure;
	long x;
	long y;
	bool kept;
	bool locked;
	long structureToBecaume;
	long newIndex;
} wsFusion2D;


typedef struct
{
	double val;
	unsigned long x;
	unsigned long y;
} maxStruct2D;

typedef struct
{
	double xx;
	double yy;
	double xy;
} varianceStruct2D;

typedef struct
{
	double real;
	double imag;
} cpx;

typedef struct
{
	double max;
	double mean;
	double angleMax;
} gaborCoefs;

typedef struct
{
	unsigned long peri;
	unsigned long surf;
	maxStruct2D max;
	maxStruct2D min;
	double stdev;
	double morpho;
	val2D moment;
	val2D momentSurf;
	double mean;
	double intensity;
	varianceStruct2D var;
	val2D sigma;
	double angle;
	unsigned long index;
	unsigned long time;
	long maxLength;
	unsigned long zOrderComputed;
	cpx zernikeMoments[36];
	gaborCoefs gaborMoments[MAX_NB_FREQ_GABOR];
	bool isPseudoZernike;
	coor2D centerZernike;
	box2D boundingBox;
	unsigned short nbSubstructures;
} InfoDataRich2D;


typedef struct
{
	unsigned long peri;
	unsigned long surf;
	maxStruct2D max;
	double morpho;
	val2D moment;
	val2D momentSurf;
	double mean;
	double intensity;
	varianceStruct2D var;
	val2D sigma;
	double angle;
	unsigned long index;
	unsigned long time;
	long maxLength;
} InfoData2D;


class InfoStruct2D
{
public:
	InfoStruct2D(int index);
	~InfoStruct2D();
	InfoStruct2D* next;

	InfoDataRich2D data;
};


#endif
