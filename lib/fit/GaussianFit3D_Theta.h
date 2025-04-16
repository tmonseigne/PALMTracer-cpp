#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit3D_Theta final : GaussianFit
{
public:
	static constexpr size_t N_PARAM = 7; // I0, x0, y0, sigmaX, sigmaY, theta, offset
	static constexpr size_t N_PARAM2 = N_PARAM * N_PARAM;

	std::array<double, N_PARAM> pMin_, p_, pMax_; // Possible maintenant

	GaussianFit3D_Theta()
		: pMin_(), p_(), pMax_()
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(const unsigned short* image, size_t sizeX, size_t sizeY, double& x, double& y, double& intensity,
	         double& sigmaX, double& sigmaY, double& theta, double& offset, double& rss, size_t iterNum, size_t roi);

protected:
	void initP(double x, double y, double sigmaX, double sigmaY, double theta);
	double lm_matx(const std::vector<size_t>& t, const double* p, const std::vector<double>& dat, const std::vector<double>& weight,
	               double* alpha, double* beta, std::vector<double>& hat) override;
};
