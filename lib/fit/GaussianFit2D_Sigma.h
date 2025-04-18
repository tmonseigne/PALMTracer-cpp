#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit2D_Sigma final : GaussianFit
{
public :
	static constexpr size_t N_PARAM = 5; // I0, x0, y0, sigma, offset
	static constexpr size_t N_PARAM2 = N_PARAM * N_PARAM;

	std::array<double, N_PARAM> pMin_, p_, pMax_; // Possible maintenant

	GaussianFit2D_Sigma()
		: pMin_(), p_(), pMax_()
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(const unsigned short* image, size_t sizeX, size_t sizeY, double& x, double& y, double& intensity,
	         double& sigma, double& offset, double& rss, size_t roi);

protected:
	double lm_matx(const std::vector<size_t>& t, const double* p, const std::vector<double>& dat, const std::vector<double>& weight,
	               double* alpha, double* beta, std::vector<double>& hat) override;
	void initP(double x, double y, double sigma);
};
