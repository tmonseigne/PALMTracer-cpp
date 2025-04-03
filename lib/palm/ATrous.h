#pragma once
#include <memory.h>
#include "tools/allDataLib.h"


class ATrous
{
	FltArray temp2_;

public:
	FltArray wDecomp_;
	ATrous() = default;
	~ATrous() = default;

	//2D transform:
private:
	static int offsetTest(int i, int j, int maxSizeI, int maxSizeJ);
	static int offset(int i, int j, int maxSizeJ);
	static int testSizes2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	static void projectionV02D(unsigned short* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY);
	int b3SplineFiltering2D(float* srcData, float* destData, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long currentPlane);

public:
	int ATrousDecomposition2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
};

