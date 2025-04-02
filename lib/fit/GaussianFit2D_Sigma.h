#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit2D_Sigma : GaussianFit
{
public :
	static constexpr unsigned short N_PARAM = 5;
	static constexpr unsigned short N_PARAM2 = N_PARAM * N_PARAM;

	std::array<float, N_PARAM> pMin_, p_, pMax_;  // Possible maintenant

	GaussianFit2D_Sigma() { nParam_ = N_PARAM; nParam2_ = N_PARAM2; } // I0, x0, y0, sigma, offset

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS,
	         int iterNum, unsigned short size);

protected:
	float lm_matx(std::vector<unsigned int>& t/*IN*/, float* p/*IN*/, std::vector<float>& y_dat/*IN*/, std::vector<float>& weight/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/,
	              std::vector<float>& y_hat/*OUT*/) override;
	void initP(float x, float y, float sigma);
};
