#pragma once
#include "GaussianFit.h"


class GaussianFit3D : GaussianFit
{
protected:
	float deriv1I0(float x, float y, float x0, float y0, float sigmaX, float sigmaY, float offset) const;
	float deriv1x0(float x, float y, float i0, float x0, float y0, float sigmaX, float sigmaY, float offset) const;
	float deriv1y0(float x, float y, float i0, float x0, float y0, float sigmaX, float sigmaY, float offset) const;
	float deriv1sigmaX(float x, float y, float i0, float x0, float y0, float sigmaX, float sigmaY, float offset) const;
	float deriv1sigmaY(float x, float y, float i0, float x0, float y0, float sigmaX, float sigmaY, float offset) const;
	float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) override;
	void inverse(float* MM, float* Minv) override;

public:
	float theta_;

	GaussianFit3D() : theta_(0) {}

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& sigmaY, float& angle, float& offset, float& RSS, int iterNum, unsigned short size);
};
