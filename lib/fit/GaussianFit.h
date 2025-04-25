#pragma once
#include <vector>

enum class fit_mode: uint8_t { gauss_2d, gauss_2d_sigma, gauss_3d, gauss_3d_theta };


class GaussianFit
{
public:
	explicit GaussianFit(fit_mode mode = fit_mode::gauss_2d);
	~GaussianFit() = default;
	void Run(const unsigned short* image, size_t sizeX, size_t sizeY, double& x, double& y, double& intensity,
	         double& sigmaX, double& sigmaY, double& theta, double& offset, double& rss, size_t roi);

protected:
	fit_mode mode_;
	size_t n_, nParams_, nParams2_;
	double sigma_, theta_;
	std::vector<double> pMin_, p_, pMax_;
	static constexpr size_t N_ITER = 50; // For X,Y,Theta le nombre d'iteration était de 140 et sinon à 10.
	static constexpr double EPSILON_1 = 1e-6, EPSILON_2 = 1e-6, EPSILON_3 = 1e-5, EPSILON_4 = 1e-2;
	static constexpr double LAMBDA_0 = 1e-2, LAMBDA_EPSILON_MAX = 1e7, LAMBDA_EPSILON_MIN = 1e-7, LAMBDA_UP = 11, LAMBDA_DN = 9;

	//==================== Run Computes ====================
	static void preComputeFit(double& min, double& max, const unsigned short* image, size_t width, size_t x, size_t y,
	                          std::vector<double>& dat, std::vector<double>& weight, std::vector<size_t>& t, size_t roi);
	void initP(double x, double y, double sigmaX, double sigmaY, double theta, double minI, double maxI);
	double computeRSS(std::vector<double>& dat, std::vector<double>& hat) const;

	//==================== PSF ====================
	double psf(const std::vector<double>& p, double x, double y) const;
	static double psf(double x, double y, double i0, double sigma, double offset);
	static double psf(double x, double y, double i0, double sigmaX, double sigmaY, double theta, double offset);

	//==================== Matrix ====================
	void deriv(double* out, const std::vector<double>& p, double x, double y) const;
	static void deriv(double* out, double x, double y, double i0, double sigma, bool has_sigma);
	static void deriv(double* out, double x, double y, double i0, double sigmaX, double sigmaY, double theta, bool has_theta);
	double chi2(std::vector<double>& alpha, std::vector<double>& beta, const std::vector<double>& dydp,
	            const std::vector<double>& weight, const std::vector<double>& deltaY) const;
	double lm_matx(const std::vector<size_t>& t, const std::vector<double>& p, const std::vector<double>& dat, const std::vector<double>& weight,
	               std::vector<double>& alpha, std::vector<double>& beta, std::vector<double>& hat) const;
	void inverse(std::vector<double>& in, std::vector<double>& out) const;
};
