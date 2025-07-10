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

//-------------------------------------------------------------------------------------------------
void Deriv(double* out, const double i0, const double x, const double y, const double sigma, const bool has_sigma)
{
	const auto X_2 = x * x, Y_2 = y * y, Sigma_2 = sigma * sigma;
	const auto Exp = exp(-(X_2 + Y_2) / (2 * Sigma_2));
	const auto P = i0 * Exp;

	out[0] = 1;                                                      // offset	: ∂PSF/∂offset
	out[1] = Exp;                                                    // I0		: ∂PSF/∂I0
	out[2] = P * x / Sigma_2;                                        // x0		: ∂PSF/∂x0
	out[3] = P * y / Sigma_2;                                        // y0		: ∂PSF/∂y0
	if (has_sigma) { out[4] = P / (sigma * Sigma_2) * (X_2 + Y_2); } // sigma	: ∂PSF/∂sigma
}

//-------------------------------------------------------------------------------------------------
void Deriv(double* out, const double i0, const double x, const double y, const double sigmaX, const double sigmaY,
	const double theta, const bool has_theta)
{
	// Calculs préliminaires
	const auto X_2 = x * x, Y_2 = y * y, X_Y = x * y, SigmaX_2 = sigmaX * sigmaX, SigmaY_2 = sigmaY * sigmaY;

	const auto CosT = cos(theta), SinT = sin(theta);
	const auto CosT_2 = CosT * CosT, SinT_2 = SinT * SinT, CosSinT = CosT * SinT,
		XCosT = x * CosT, XSinT = x * SinT, YCosT = y * CosT, YSinT = y * SinT;

	const auto Factor1 = (XCosT - YSinT) / sigmaX,
		Factor2 = (XSinT + YCosT) / sigmaY,
		Factor3 = 2 * X_Y * CosSinT,
		Factor4 = 2 * X_Y / SigmaX_2 - 2 * X_Y / SigmaY_2;
	const auto Exp = exp(-(Factor1 * Factor1 + Factor2 * Factor2));
	const auto P = i0 * Exp;

	out[0] = 1;                                                                                                                              // offset
	out[1] = Exp;                                                                                                                            // I0
	out[2] = P * (-((-2 * CosT_2 * x + 2 * y * CosSinT) / SigmaX_2) - (-2 * SinT_2 * x - 2 * y * CosSinT) / SigmaY_2);                       // x0
	out[3] = P * (-((2 * SinT_2 * y - 2 * x * CosSinT) / SigmaX_2) - (-2 * CosT_2 * y - 2 * x * CosSinT) / SigmaY_2);                        // y0
	out[4] = P * (2 * (XCosT * XCosT - Factor3 + YSinT * YSinT) / (sigmaX * SigmaX_2));                                                      // sigma X
	out[5] = P * (2 * (XSinT * XSinT + Factor3 + YCosT * YCosT) / (sigmaY * SigmaY_2));                                                      // sigma Y
	if (has_theta) { out[6] = P * (2 * CosSinT * ((X_2 - Y_2) / SigmaX_2 + (Y_2 - X_2) / SigmaY_2) - SinT_2 * Factor4 + CosT_2 * Factor4); } // theta
}

} // namespace palm_tracer
