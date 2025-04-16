/*******************************************************************************************************************
								Link class to treat data from segmentation of images 3D (16bit only)
*******************************************************************************************************************/
#pragma once


using val3D = struct
{
	double x, y, z;
};

using coor3D = struct
{
	long x;
	long y;
	long z;
};

using wsFusion3D = struct
{
	float maxPotentialWithNeighbor;
	long structureInContact;
	float minPotentialOnTheStructure;
	float maxPotentialOnTheStructure;
	long x;
	long y;
	long z;
	bool kept;
	bool locked;
	long structureToBecaume;
	long newIndex;
};


using anglesEuler = struct
{
	double phi;
	double teta;
	double psi;
};


using maxStruct3D = struct
{
	double val;
	unsigned long x;
	unsigned long y;
	unsigned long z;
};

using varianceStruct3D = struct
{
	double xx;
	double yy;
	double zz;
	double xy;
	double yz;
	double xz;
};


using InfoData3D = struct
{
	unsigned long surfVol;
	unsigned long volume;
	maxStruct3D max;
	val3D moment;
	val3D momentSurf;
	double mean;
	double intensity;
	varianceStruct3D var;
	val3D sigma;
	anglesEuler angle;
	double det;
	unsigned long index;
	unsigned long time;
};


class InfoStruct3D
{
public:
	InfoStruct3D(const int index)
	{
		data.index = index;
		next = nullptr;
	}

	~InfoStruct3D() = default;
	InfoStruct3D* next;

	InfoData3D data;
};
