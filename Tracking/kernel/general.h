/*******************************************************************************************************************
								General Structures used for the Tracking
*******************************************************************************************************************/

#ifndef __general__
#define __general__
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>
#include <math.h>
#include <wchar.h>

//#include "print_debug.h"
#include "../HTML/print.h"
#define	MIN(i,j)	((i)<(j) ? (i) : (j) )


//void getOutputName(char* in, char* out);
int getFileName(wchar_t* folder, wchar_t* input, wchar_t* output);

class PublicParameters
{
public:
	PublicParameters() { };
	~PublicParameters() { };
	double distMaxNeighbour;
	double distMaxPastAndFuture;
	double dZDX;
	int nbNodeAllocated;
	int nbNodeDeallocated;
	int nbChainNeighborAllocated;
	int nbChainNeighborDeallocated;
	wchar_t statFile[1024];
	wchar_t** segImagesIn;
	long segImagesInCount;
	wchar_t fileRes[1024];
	wchar_t folderOut[1024];
	wchar_t fileResASCII[1024];
	wchar_t fileArg0[1024];
	long model; //the number of the used model
	double likenessIntensityMin;
	double temperatureDecrease;
	long nbDim;
	long nbMinSuccessiveStructures;
	bool allowFusionAndFission;
	long nbTemporalWindows;
	double costBirthDeath;

	bool drawSpeedVector2D;

	double stdIntensity;
	double stdMean;
	double stdDistance;
	double covMeanIntensity;

	double maxX;
	double maxY;
	double maxZ;

	// new members
	double* pointsList;
	double* trc;
	unsigned long pointsListSize;
	unsigned long trcSize;
	unsigned long nbTimes;
	long nbTraj;
	int NbSegParams;

	int NbIterations;
};


#endif	//__general__
