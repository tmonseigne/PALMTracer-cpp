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

	int ATrousDecomposition2D(unsigned short* image, size_t height, size_t width, bool* planesToUse, size_t nbPlane, bool realProj);

	//2D transform:
private:
	static int offsetTest(int i, int j, int maxSizeI, int maxSizeJ);
	static int offset(int i, int j, int maxSizeJ);
	int testSizes2D(size_t height, size_t width, bool* planesToUse, size_t& nbPlane);
	void projectionV02D(const unsigned short* image, size_t realProj, size_t height, size_t width);
	int b3SplineFiltering2D(const float* srcData, float* destData, size_t height, size_t width, size_t currentPlane);
};

