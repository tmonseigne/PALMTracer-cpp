#pragma once
#include <vector>

class GaussianFit
{
public:
	size_t n_, nParam_, nParam2_;
	static constexpr size_t N_ITER = 10;
	static constexpr double EPSILON_1 = 1e-6, EPSILON_2 = 1e-6, EPSILON_3 = 1e-5, EPSILON_4 = 1e-2;
	static constexpr double LAMBDA_0 = 1e-2, LAMBDA_EPSILON_MAX = 1e7, LAMBDA_EPSILON_MIN = 1e-7, LAMBDA_UP = 11, LAMBDA_DN = 9;

	GaussianFit() : n_(0), nParam_(0), nParam2_(0) { }
	virtual ~GaussianFit() = default;

protected:
	static double psf(double x, double y, double i0, double sigma, double offset);
	static double psf(double x, double y, double i0, double sigmaX, double sigmaY, double theta, double offset);

	static void deriv(double* out, double x, double y, double i0, double sigma, bool has_sigma);
	static void deriv(double* out, double x, double y, double i0, double sigmaX, double sigmaY, double theta, bool has_theta);

	double chi2(double* alpha, double* beta, const std::vector<double>& dydp, const std::vector<double>& weight, const std::vector<double>& deltaY) const;
	double computeRSS(std::vector<double>& dat, std::vector<double>& hat) const;

	virtual double lm_matx(const std::vector<size_t>& t, const double* p, const std::vector<double>& dat, const std::vector<double>& weight,
	                       double* alpha, double* beta, std::vector<double>& hat) = 0;

	void inverse(double* in, double* out) const;

	static void preComputeFit(double& min, double& max, const unsigned short* image, size_t width, size_t x, size_t y,
	                          std::vector<double>& dat, std::vector<double>& weight, std::vector<size_t>& t, size_t roi);
};
