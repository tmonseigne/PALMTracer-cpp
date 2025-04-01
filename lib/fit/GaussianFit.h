#pragma once

class GaussianFit
{
public:
	virtual ~GaussianFit() = default;

private:
	virtual float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) = 0;
	virtual void inverse(float* MM, float* Minv) { }
};
