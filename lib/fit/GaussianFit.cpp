#include "GaussianFit.h"

#include <algorithm>
#include <cmath>
#include <execution>
#include <numeric>
#include <vector>

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

void GaussianFit::deriv(float* out, const float x, const float y, const float i0, const float sigma, const bool has_sigma)
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

void GaussianFit::deriv(float* out, const float x, const float y, const float i0, const float sigmaX, const float sigmaY, const float theta, const bool has_theta)
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

float GaussianFit::chi2(float* alpha, float* beta, const std::vector<float>& dydp, const std::vector<float>& weight, const std::vector<float>& deltaY) const
{
	std::vector<float> tmp(static_cast<size_t>(n_ * nParam_));
	for (int i = 0; i < n_; ++i) { for (int j = 0; j < nParam_; ++j) { tmp[i * nParam_ + j] = weight[i] * dydp[i * nParam_ + j]; } }

	// Computing alpha (A v�rifier)
	for (int i = 0; i < nParam_; ++i) {
		for (int j = 0; j < nParam_; ++j) {
			alpha[nParam_ * i + j] = 0;
			for (int k = 0; k < n_; k++) { alpha[nParam_ * i + j] += dydp[k * nParam_ + i] * tmp[k * nParam_ + j]; }
		}
	}

	// Computing beta
	tmp.resize(n_);
	std::ranges::transform(weight, deltaY, tmp.begin(), std::multiplies<>{});

	for (int i = 0; i < nParam_; ++i) {
		beta[i] = 0;
		for (int j = 0; j < n_; ++j) { beta[i] += dydp[j * nParam_ + i] * tmp[j]; }
	}

	// Computing Chi_sq
	return std::transform_reduce(std::execution::seq, deltaY.begin(), deltaY.end(), tmp.begin(), 0.0f);
}

float GaussianFit::computeRSS(std::vector<float>& dat, std::vector<float>& hat) const
{
	// computing "dat" mean & "hat" mean
	const float Dat_Mean = std::accumulate(dat.begin(), dat.end(), 0.0f) / static_cast<float>(n_),
	            Hat_Mean = std::accumulate(hat.begin(), hat.end(), 0.0f) / static_cast<float>(n_);

	// Normalisation en place
	for (int i = 0; i < n_; ++i) {
		dat[i] -= Dat_Mean;
		hat[i] -= Hat_Mean;
	}

	// Calcul des t1, t2, t3
	const float T1 = std::inner_product(dat.begin(), dat.end(), hat.begin(), 0.0f), // dat * hat
	            T2 = std::inner_product(dat.begin(), dat.end(), dat.begin(), 0.0f), // dat * dat
	            T3 = std::inner_product(hat.begin(), hat.end(), hat.begin(), 0.0f), // hat * hat
	            R = T1 / std::sqrt(T2 * T3);                                        // Calcul du R

	return std::isnan(R) ? 0 : R * R;
}

void GaussianFit::inverse(float* in, float* out) const
{
	std::vector<float> m1(in, in + nParam2_); // Copie de in dans M1
	std::vector<float> mi1(nParam2_, 0.0f);
	for (int i = 0; i < nParam_; ++i) { mi1[i * nParam_ + i] = 1.0f; } // Initialisation de la matrice identit�

	// �limination de Gauss-Jordan
	for (int i = 0; i < nParam_; ++i) {
		const float Pivot = m1[i * nParam_ + i];
		const float Pivot_Inv = 1.0f / Pivot;

		// Normalisation de la ligne i
		for (int j = 0; j < nParam_; ++j) {
			m1[i * nParam_ + j] *= Pivot_Inv;
			mi1[i * nParam_ + j] *= Pivot_Inv;
		}

		// Mise � z�ro des autres lignes
		for (int j = 0; j < nParam_; ++j) {
			if (j != i) {
				const float Factor = m1[j * nParam_ + i];

				for (int k = 0; k < nParam_; ++k) {
					m1[j * nParam_ + k] -= Factor * m1[i * nParam_ + k];
					mi1[j * nParam_ + k] -= Factor * mi1[i * nParam_ + k];
				}
			}
		}
	}

	// Copie finale de Mi1 vers out
	std::ranges::copy(mi1, out);
}

void GaussianFit::preComputeFit(float& min, float& max, const unsigned short* image, const unsigned long width,
                                const unsigned int x, const unsigned int y, std::vector<float>& dat, std::vector<float>& weight, std::vector<unsigned int>& t,
                                const unsigned short roi)
{
	min = std::numeric_limits<float>::max();    // Plus grand float possible (FLT_MAX) // Not present in GPU original file why ?
	max = std::numeric_limits<float>::lowest(); // Plus petit float possible (-FLT_MAX) // Not present in GPU original file why ?

	const unsigned short Half = roi / 2;
	for (int i = -Half; i <= Half; ++i) {
		const int Idx = Half + (i + Half) * roi;
		for (int j = -Half; j <= Half; ++j) {
			dat[j + Idx] = static_cast<float>(image[x + j + (y + i) * width]);
			const auto TIdx = (j + Idx) * 2;
			t[TIdx] = x + j;
			t[TIdx + 1] = y + i;

			max = std::max(max, dat[j + Idx]); // Not present in GPU original file why ?
			min = std::min(min, dat[j + Idx]); // Not present in GPU original file why ?
		}
	}

	std::ranges::fill(weight, 1);
	//for (int i = 0; i < n_; ++i) {
	//	weight[i] = 1; //sqrt((n_-nParam_+1) / (dat[ii]*dat[ii]))
	//	//weight[i] *= weight[i];
	//}
}
