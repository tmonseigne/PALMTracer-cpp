#pragma once
#include <vector>

constexpr float EPSILON_1 = 1e-6f;
constexpr float EPSILON_2 = 1e-6f;
constexpr float EPSILON_3 = 1e-5f;
constexpr float EPSILON_4 = 1e-2f;
constexpr float LAMBDA_0 = 1e-2f;
constexpr float LAMBDA_UP_FAC = 11;
constexpr float LAMBDA_DN_FAC = 9;

class GaussianFit
{
public:
	unsigned short n_;
	unsigned short nParam_;
	std::vector<float> pMin_, p_, pMax_;

	GaussianFit() : n_(0), nParam_(0) { }
	virtual ~GaussianFit() = default;

protected:
	static float psf(float x, float y, float i0, float sigma, float offset);
	static float psf(float x, float y, float i0, float sigmaX, float sigmaY, float theta, float offset);

	static void deriv(float* out, float x, float y, float i0, float sigma, bool has_sigma);
	static void deriv(float* out, float x, float y, float i0, float sigmaX, float sigmaY, float theta, bool has_theta);

	float chi2(float* alpha, float* beta, const float* dydp, const float* weight, const float* deltaY) const;

	virtual float lm_matx(unsigned int* t, float* p, float* y_dat, float* weight, float* alpha, float* beta, float* y_hat) = 0;
	void inverse(float* in, float* out) const;

	void preComputeFit(float& min, float& max, const unsigned short* image, unsigned long width, unsigned int x, unsigned int y,
	                   float* y_dat, float* weight, unsigned int* t, unsigned short roi) const;
};
