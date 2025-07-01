#include "PSF.hpp"

#include <cmath>
#include <execution>

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
double PSF(const double offset, const double i0, const double x, const double y, const double sigma)
{
	return i0 * exp(-(x * x + y * y) / (2 * (sigma * sigma))) + offset;
}

//-------------------------------------------------------------------------------------------------
double PSF(const double offset, const double i0, const double x, const double y, const double sigmaX, const double sigmaY, const double theta)
{
	const auto CosTheta = cos(theta), SinTheta = sin(theta);
	const auto Factor1  = (x * CosTheta - y * SinTheta) / sigmaX,
	           Factor2  = (x * SinTheta + y * CosTheta) / sigmaY;

	return i0 * exp(-(Factor1 * Factor1 + Factor2 * Factor2)) + offset;
}
} // namespace palm_tracer
