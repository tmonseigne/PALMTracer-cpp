#pragma once
#include <vector>

constexpr double EPSILON_1 = 1e-6;
constexpr double EPSILON_2 = 1e-6;
constexpr double EPSILON_3 = 1e-5;
constexpr double EPSILON_4 = 1e-2;
constexpr double LAMBDA_0 = 1e-2;
constexpr double LAMBDA_EPSILON_MAX = 1e7;
constexpr double LAMBDA_EPSILON_MIN = 1e-7;
constexpr double LAMBDA_UP_FAC = 11;
constexpr double LAMBDA_DN_FAC = 9;

class GaussianFit
{
public:
	size_t n_, nParam_, nParam2_;

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
