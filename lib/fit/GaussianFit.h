#pragma once

class GaussianFit
{
public:
	virtual ~GaussianFit() = default;

protected:
	static float psf(float x, float y, float i0, float sigma, float offset);
	static float psf(float x, float y, float i0, float sigmaX, float sigmaY, float theta, float offset);
	virtual float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) = 0;
	virtual void inverse(float* MM, float* Minv) { }
};
