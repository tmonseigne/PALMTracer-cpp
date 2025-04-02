#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit3D_Theta final : GaussianFit
{
public:
	static constexpr unsigned short N_PARAM = 7; // I0, x0, y0, sigmaX, sigmaY, theta, offset
	static constexpr unsigned short N_PARAM2 = N_PARAM * N_PARAM;

	std::array<float, N_PARAM> pMin_, p_, pMax_; // Possible maintenant

	GaussianFit3D_Theta(): pMin_(), p_(), pMax_()
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& i0,
	         float& sigmaX, float& sigmaY, float& angle, float& offset, float& rss, int iterNum, unsigned short size);

protected:
	void initP(float x, float y, float sigmaX, float sigmaY, float theta);
	float lm_matx(const std::vector<unsigned int>& t, const float* p, const std::vector<float>& dat, const std::vector<float>& weight,
	              float* alpha, float* beta, std::vector<float>& hat) override;
};
