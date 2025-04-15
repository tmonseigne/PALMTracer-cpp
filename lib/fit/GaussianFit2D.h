#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit2D final : GaussianFit
{
public:
	static constexpr size_t N_PARAM = 4; // I0, x0, y0, offset
	static constexpr size_t N_PARAM2 = N_PARAM * N_PARAM;

	std::array<double, N_PARAM> pMin_, p_, pMax_; // Possible maintenant

	double sigma_;

	GaussianFit2D() : pMin_(), p_(), pMax_(), sigma_(0)
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(const unsigned short* image, size_t sizeX, size_t sizeY, double& x, double& y, double& intensity,
	         const double& sigma, double& offset, double& rss, size_t iterNum, size_t roi);

protected:
	void initP(double x, double y);
	double lm_matx(const std::vector<size_t>& t, const double* p, const std::vector<double>& dat, const std::vector<double>& weight,
	              double* alpha, double* beta, std::vector<double>& hat) override;
};
