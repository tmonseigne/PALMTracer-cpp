#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit3D final : GaussianFit
{
public:
	static constexpr size_t N_PARAM = 6; // I0, x0, y0, sigmaX, sigmaY, offset
	static constexpr size_t N_PARAM2 = N_PARAM * N_PARAM;

	std::array<double, N_PARAM> pMin_, p_, pMax_; // Possible maintenant

	double theta_;

	GaussianFit3D()
		: pMin_(), p_(), pMax_(), theta_(0)
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(const unsigned short* image, size_t sizeX, size_t sizeY, double& x, double& y, double& intensity,
	         double& sigmaX, double& sigmaY, const double& theta, double& offset, double& rss, size_t roi);

protected:
	void initP(double x, double y, double sigmaX, double sigmaY);
	double lm_matx(const std::vector<size_t>& t, const double* p, const std::vector<double>& dat, const std::vector<double>& weight,
	               double* alpha, double* beta, std::vector<double>& hat) override;
};
