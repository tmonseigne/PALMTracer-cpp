#include "GaussianFit.h"

#include <cmath>

float GaussianFit::psf(const float x, const float y, const float i0, const float sigma, const float offset)
{
	return i0 * exp(-(x * x + y * y) / (2 * (sigma * sigma))) + offset;
}

float GaussianFit::psf(const float x, const float y, const float i0, const float sigmaX, const float sigmaY, const float theta, const float offset)
{
	const auto CosTheta = cos(theta), SinTheta = sin(theta);
	const auto Factor1 = (x * CosTheta - y * SinTheta) / sigmaX,
	           Factor2 = (x * SinTheta + y * CosTheta) / sigmaY;

	return i0 * exp(-(Factor1 * Factor1 + Factor2 * Factor2)) + offset;
}
