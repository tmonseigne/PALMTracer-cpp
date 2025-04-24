/*******************************************************************************************************************
								Link class to treat data from segmentation of images 3D (16bit only)
*******************************************************************************************************************/
#pragma once


using val_3d = struct
{
	double x, y, z;
};

using coor_3d = struct
{
	int x, y, z;
};

using ws_fusion_3d = struct
{
	float maxPotentialWithNeighbor;
	int structureInContact;
	float minPotentialOnTheStructure;
	float maxPotentialOnTheStructure;
	int x, y, z;
	bool kept;
	bool locked;
	int structToBecome;
	int newIndex;
};


using angles_euler = struct
{
	double phi, teta, psi;
};


using max_3d = struct
{
	double val;
	size_t x, y, z;
};

using variance_3d = struct
{
	double xx, yy, zz, xy, yz, xz;
};


using info_3d = struct
{
	size_t surfVol;
	size_t volume;
	max_3d max;
	val_3d moment;
	val_3d momentSurf;
	double mean;
	double intensity;
	variance_3d var;
	val_3d sigma;
	angles_euler angle;
	double det;
	size_t index;
	size_t time;
};


class Info3D
{
public:
	Info3D(const int index)
	{
		data_.index = index;
		next_ = nullptr;
	}

	~Info3D() = default;
	Info3D* next_;

	info_3d data_;
};
