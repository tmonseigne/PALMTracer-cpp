#pragma once
#include <array>
#include "GaussianFit.h"


class GaussianFit3D : GaussianFit
{
public:
	static constexpr unsigned short N_PARAM = 6; // I0, x0, y0, sigmaX, sigmaY, offset
	static constexpr unsigned short N_PARAM2 = N_PARAM * N_PARAM;

	std::array<float, N_PARAM> pMin_, p_, pMax_; // Possible maintenant

	float theta_;

	GaussianFit3D() : pMin_(), p_(), pMax_(), theta_(0)
	{
		nParam_ = N_PARAM;
		nParam2_ = N_PARAM2;
	}

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& sigmaY, float& angle,
	         float& offset, float& RSS, int iterNum, unsigned short size);

protected:
	void initP(float x, float y, float sigmaX, float sigmaY);
	float lm_matx(const std::vector<unsigned int>& t, const float* p, const std::vector<float>& dat, const std::vector<float>& weight, float* alpha, float* beta, std::vector<float>& hat) override;
};
