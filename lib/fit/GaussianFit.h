#pragma once

class GaussianFit
{
public:
	unsigned short dimFit_;
	unsigned short n_;
	unsigned short nParam_;
	float ref_;

	GaussianFit() : dimFit_(0), n_(0), nParam_(0), ref_(0) { }
	virtual ~GaussianFit() = default;

protected:
	static float psf(float x, float y, float i0, float sigma, float offset);
	static float psf(float x, float y, float i0, float sigmaX, float sigmaY, float theta, float offset);

	static void deriv(float* out, float x, float y, float i0, float sigma, bool has_sigma);
	static void deriv(float* out, float x, float y, float i0, float sigmaX, float sigmaY, float theta, bool has_theta);

	float chi2(float* alpha, float* beta, const float* dydp, const float* weight, const float* deltaY) const;

	void inverse(float* in, float* out) const;

	virtual float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) = 0;
};
