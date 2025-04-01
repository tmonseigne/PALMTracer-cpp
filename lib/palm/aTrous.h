#pragma once
#include <memory.h>
#include "tools/allDataLib.h"
#include "tools/miscTools.h"


class ATrous
{
	//general part:
private:
	FltArray* temp1;
	FltArray* temp2;
	FltArray* inputV0;
	FltArray* VN;
	FltArray* VNplus1;

public:
	FltArray* wDecomp;
	ATrous();
	~ATrous();

	//tools:
public:
	double stdRef(long dim, unsigned long realProj, unsigned long plan);

	//2D transform:
private:
	int offsetTest(int i, int j, int maxSizeI, int maxSizeJ);
	int offset(int i, int j, int maxSizeJ);
	int testSizes2D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	int testSizes2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	void projectionV02D(unsigned short* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY);
	void projectionV02D(float* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY);
	int b3SplineFiltering2D(float* srcData, float* destData, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long currentPlane);

public:
	int aTrousDecomposition2D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
	int aTrousDecomposition2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
};

