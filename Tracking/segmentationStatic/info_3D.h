

/*******************************************************************************************************************
								Link class to treat data from segmentation of images 3D (16bit only)
*******************************************************************************************************************/

#ifndef __link_3D__
#define __link_3D__

#include <stdio.h>
#include "tools/allDataLib.h"
extern int __nbInfoStrucAlloc;
extern int __nbInfoStrucDealloc;

typedef struct {
	double	x;
	double	y;
	double	z;
}	val3D; 

typedef struct {
	long	x;
	long	y;
	long	z;
}	coor3D;

typedef struct {
	float	maxPotentialWithNeighbor; 
	long	structureInContact;
	float	minPotentialOnTheStructure;
	float	maxPotentialOnTheStructure;
	long x;
	long y;
	long z;
	bool kept;
	bool locked;
	long structureToBecaume;
	long newIndex;
}	wsFusion3D;


typedef struct {
	double	phi;
	double	teta;
	double	psi;
}	anglesEuler;


typedef struct {
	double			val;
	unsigned long	x;
	unsigned long	y;
	unsigned long	z;
}	maxStruct3D;

typedef struct {
	double	xx;
	double	yy;
	double	zz;
	double	xy;
	double	yz;
	double	xz;
}	varianceStruct3D;


typedef struct {
	unsigned long		surfVol;
	unsigned long		volume;
	maxStruct3D			max;
	double				morpho;
	val3D				moment;
	val3D				momentSurf;
	double				mean;
	double				intensity;
	varianceStruct3D	var;
	val3D				sigma;
	anglesEuler			angle;
	double				det;
	unsigned long		index;
	unsigned long		time;
}	InfoData3D;


class InfoStruct3D{
public:
	InfoStruct3D(int index);
	~InfoStruct3D();
	InfoStruct3D* next;

	InfoData3D data;
	
};



#endif
