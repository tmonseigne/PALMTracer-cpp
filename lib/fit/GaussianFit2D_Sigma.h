#pragma once
#include "GaussianFit.h"


class GaussianFit2D_Sigma : GaussianFit
{
protected:
	static float deriv1I0(float x, float y, float x0, float y0, float sigma);
	static float deriv1x0(float x, float y, float i0, float x0, float y0, float sigma);
	static float deriv1y0(float x, float y, float i0, float x0, float y0, float sigma);
	static float deriv1sigma(float x, float y, float i0, float x0, float y0, float sigma);
	float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) override;
	void inverse(float* MM, float* Minv) override;

public :
	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS, int iterNum, unsigned short size);
};
