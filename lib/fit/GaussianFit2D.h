#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit2D : GaussianFit
{
public:
	static constexpr unsigned short N_PARAM = 4;
	static constexpr unsigned short N_PARAM2 = N_PARAM * N_PARAM;

	std::array<float, N_PARAM> pMin_, p_, pMax_;  // Possible maintenant

	float sigma_;

	GaussianFit2D() : sigma_(0) { nParam_ = N_PARAM; nParam2_ = N_PARAM2; } // I0, x0, y0, offset

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS,
	         int iterNum, unsigned short size);

protected:
	void initP(float x, float y);
	float lm_matx(std::vector<unsigned int>& t/*IN*/, float* p/*IN*/, std::vector<float>& y_dat/*IN*/, std::vector<float>& weight/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/,
	              std::vector<float>& y_hat/*OUT*/) override;
};
