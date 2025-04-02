#pragma once
#include <vector>

constexpr float EPSILON_1 = 1e-6f;
constexpr float EPSILON_2 = 1e-6f;
constexpr float EPSILON_3 = 1e-5f;
constexpr float EPSILON_4 = 1e-2f;
constexpr float LAMBDA_0 = 1e-2f;
constexpr float LAMBDA_EPSILON = 1e7f;
constexpr float LAMBDA_UP_FAC = 11;
constexpr float LAMBDA_DN_FAC = 9;

class GaussianFit
{
public:
	unsigned short n_, nParam_, nParam2_;

	GaussianFit() : n_(0), nParam_(0), nParam2_(0) { }
	virtual ~GaussianFit() = default;

protected:
	static float psf(float x, float y, float i0, float sigma, float offset);
	static float psf(float x, float y, float i0, float sigmaX, float sigmaY, float theta, float offset);

	static void deriv(float* out, float x, float y, float i0, float sigma, bool has_sigma);
	static void deriv(float* out, float x, float y, float i0, float sigmaX, float sigmaY, float theta, bool has_theta);

	float chi2(float* alpha, float* beta, const std::vector<float>& dydp, const  std::vector<float>& weight, const std::vector<float>& deltaY) const;
	float computeRSS(std::vector<float>& dat, std::vector<float>& hat) const;

	virtual float lm_matx(const std::vector<unsigned int>& t, const  float* p, const std::vector<float>& dat, const std::vector<float>& weight, float* alpha, float* beta, std::vector<float>& hat) = 0;
	void inverse(float* in, float* out) const;

	void preComputeFit(float& min, float& max, const unsigned short* image, unsigned long width, unsigned int x, unsigned int y,
	                   std::vector<float>& dat, std::vector<float>& weight, std::vector<unsigned int>& t, unsigned short roi) const;
};
