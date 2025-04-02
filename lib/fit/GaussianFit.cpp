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

void GaussianFit::deriv(float* out, float x, float y, float i0, float sigma, bool has_sigma)
{
	const auto X_2 = x * x, Y_2 = y * y, Sigma_2 = sigma * sigma;
	const auto Exp = exp(-(X_2 + Y_2) / (2 * Sigma_2));
	const auto P = i0 * Exp;

	out[0] = Exp;             // I0
	out[1] = P * x / Sigma_2; // x0
	out[2] = P * y / Sigma_2; // y0

	if (has_sigma) {
		out[3] = P / (sigma * Sigma_2) * (X_2 + Y_2); // sigma
		out[4] = 1;                                   // offset
	}
	else { out[3] = 1; } // offset (quand sigma est fixe)
}

void GaussianFit::deriv(float* out, float x, float y, float i0, float sigmaX, float sigmaY, float theta, bool has_theta)
{
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

	out[0] = Exp;                                                                                                      // I0
	out[1] = P * (-((-2 * CosT_2 * x + 2 * y * CosSinT) / SigmaX_2) - (-2 * SinT_2 * x - 2 * y * CosSinT) / SigmaY_2); // x0
	out[2] = P * (-((2 * SinT_2 * y - 2 * x * CosSinT) / SigmaX_2) - (-2 * CosT_2 * y - 2 * x * CosSinT) / SigmaY_2);  // y0
	out[3] = P * (2 * (XCosT * XCosT - Factor3 + YSinT * YSinT) / (sigmaX * SigmaX_2));                                // sigma X
	out[4] = P * (2 * (XSinT * XSinT + Factor3 + YCosT * YCosT) / (sigmaY * SigmaY_2));                                // sigma Y

	if (has_theta) {
		out[5] = P * (2 * CosSinT * ((X_2 - Y_2) / SigmaX_2 + (Y_2 - X_2) / SigmaY_2) - SinT_2 * Factor4 + CosT_2 * Factor4); // theta
		out[6] = 1;                                                                                                           // offset
	}
	else { out[5] = 1; } // offset (quand theta est fixe)	
}

float GaussianFit::chi2(float* alpha, float* beta, const float* dydp, const float* weight, const float* deltaY, const size_t nParam)
{
	//float inter[nParam * n_];
	float* inter = (float*)malloc(nParam * n_ * sizeof(float));
	for (int elt = 0; elt < n_; elt++) { for (int i = 0; i < nParam; i++) { inter[elt * nParam + i] = weight[elt] * dydp[elt * nParam + i]; } }

	// Computing alpha (A vérifier)
	for (int j = 0; j < nParam; j++) {
		for (int i = 0; i < nParam; i++) {
			alpha[nParam * j + i] = 0;
			for (int cursor = 0; cursor < n_; cursor++) { alpha[nParam * j + i] += dydp[cursor * nParam + j] * inter[cursor * nParam + i]; }
			// dydp' * inter : dydp'[i*n_ + cursor]
		}
	}

	// Computing beta
	//float inter1[n_];
	float* inter1 = (float*)malloc(n_ * sizeof(float));

	for (int elt = 0; elt < n_; elt++) { inter1[elt] = weight[elt] * deltaY[elt]; }

	for (int j = 0; j < nParam; j++) {
		beta[j] = 0;
		for (int cursor = 0; cursor < n_; cursor++) { beta[j] += dydp[cursor * nParam + j] * inter1[cursor]; } // dydp' * inter : dydp'[j*n_ + cursor]
	}

	// Computing Chi_sq
	float res = 0.;
	for (int elt = 0; elt < n_; elt++) { res += deltaY[elt] * inter1[elt]; }

	// cleanup memory
	free(inter);
	free(inter1);

	return res;
}
