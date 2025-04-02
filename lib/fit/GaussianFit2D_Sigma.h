#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit2D_Sigma : GaussianFit
{
public :
	static constexpr unsigned short N_PARAM = 5; // I0, x0, y0, sigma, offset
	static constexpr unsigned short N_PARAM2 = N_PARAM * N_PARAM;

	std::array<float, N_PARAM> pMin_, p_, pMax_;  // Possible maintenant

	GaussianFit2D_Sigma(): pMin_(), p_(), pMax_()
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS,
	         int iterNum, unsigned short size);

protected:
	float lm_matx(const std::vector<unsigned int>& t, const float* p, const std::vector<float>& dat, const std::vector<float>& weight, float* alpha, float* beta,
	              std::vector<float>& hat) override;
	void initP(float x, float y, float sigma);
};
