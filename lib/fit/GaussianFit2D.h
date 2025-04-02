#pragma once
#include "GaussianFit.h"


class GaussianFit2D : GaussianFit
{
protected:
	float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) override;

public :
	float sigma_;

	GaussianFit2D() : sigma_(0) { nParam_ = 4; } // I0, x0, y0, offset

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS, int iterNum, unsigned short size);
};
