#pragma once
#include <vector>

#include "tools/allDataLib.h"


class ATrous
{
public:
	std::vector<double> tmp_;
	FltArray wDecomp_;
	ATrous() = default;
	~ATrous() = default;

	int ATrousDecomposition2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, size_t nbPlane, bool realProj);

	//2D transform:
private:
	static int offsetTest(int i, int j, int maxSizeI, int maxSizeJ);
	static int offset(int i, int j, int maxSizeJ);
	static int testSizes2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, size_t& nbPlane);
	static void projectionV02D(unsigned short* srcData, float* destData, size_t realProj, size_t height, size_t width);
	int b3SplineFiltering2D(float* srcData, float* destData, size_t height, size_t width, size_t currentPlane);
};

