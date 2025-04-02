#pragma once
#include "GaussianFit.h"


class GaussianFit2D_Sigma : GaussianFit
{
public :
	GaussianFit2D_Sigma() { nParam_ = 5; } // I0, x0, y0, sigma, offset

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS,
	         int iterNum, unsigned short size);

protected:
	float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/,
	              float* y_hat/*OUT*/) override;
	void initP(float x, float y, float sigma);
};
