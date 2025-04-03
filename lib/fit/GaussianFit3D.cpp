#include "GaussianFit3D.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>

float GaussianFit3D::lm_matx(const std::vector<unsigned int>& t, const float* p, const std::vector<float>& dat, const std::vector<float>& weight,
                             float* alpha, float* beta, std::vector<float>& hat)
{
	std::vector<float> delta_y(n_);
	for (size_t i = 0; i < n_; ++i) {
		hat[i] = psf(static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], p[3], p[4], theta_, p[5]);
		delta_y[i] = dat[i] - hat[i];
	}

	// Computing Jacobian matrix (dydp)
	std::vector<float> dydp(n_ * nParam_);
	for (size_t i = 0; i < n_; ++i) {
		deriv(&dydp[i * nParam_], static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], p[3], p[4], theta_, false);
	}

	return chi2(alpha, beta, dydp, weight, delta_y);
}

void GaussianFit3D::initP(const float x, const float y, const float sigmaX, const float sigmaY)
{
	pMin_ = {0, x - 1, y - 1, 0.1f, 0.1f, 0};
	p_ = {2000, x, y, std::numbers::sqrt2_v<float> * sigmaX, std::numbers::sqrt2_v<float> * sigmaY, 700};
	pMax_ = {65535, x + 1, y + 1, 20, 20, 65535};
}

void GaussianFit3D::Run(const unsigned short* image, const unsigned long sizeX, const unsigned long sizeY,
                        float& x, float& y, float& intensity, float& sigmaX, float& sigmaY, float& theta, float& offset,
                        float& rss, const int iterNum, unsigned short roi)
{
	if (roi % 2 == 0) { roi--; } // Not present in GPU original file why ?
	n_ = roi * roi;
	theta_ = theta;

	const unsigned short Half = roi / 2, Ux = static_cast<unsigned short>(x), Uy = static_cast<unsigned short>(y);
	if (Half <= Ux && Ux <= sizeX - Half - 1 && Half <= Uy && Uy <= sizeY - Half - 1) {
		initP(x, y, sigmaX, sigmaY);

		// Création des variables
		float lambda = LAMBDA_0;

		std::array<float, N_PARAM2> alpha, alphaR, alphaRi;
		std::array<float, N_PARAM> beta, beta_l_deltaP, p_try, deltaP;
		std::vector<float> dat(n_), hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<unsigned int> t(2 * n_);

		// Initialization of dat & t (A vérifier)
		float minI, maxI; // Not present in GPU original file why ?
		preComputeFit(minI, maxI, image, sizeX, static_cast<int>(x), static_cast<int>(y), dat, weight, t, roi);
		p_[0] = maxI - minI; // Not present in GPU original file why ?
		p_[5] = minI;        // Not present in GPU original file why ?

		float x2 = lm_matx(t, p_.data(), dat, weight, alpha.data(), beta.data(), hat);
		float x2_old = x2;

		// Trouver la valeur absolue maximale
		float m = std::ranges::max(std::views::transform(beta, [](const float b) { return std::abs(b); }));

		bool stop = false;
		if (m < EPSILON_1) { stop = true; }

		unsigned short iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			deltaP.fill(0.0f);

			std::copy_n(alpha.begin(), nParam2_, alphaR.begin()); // Copy du alpha
			for (unsigned short i = 0; i < nParam_; ++i) { alphaR[i * nParam_ + i] += lambda * alpha[i * nParam_ + i]; }
			inverse(alphaR.data(), alphaRi.data());

			for (unsigned short i = 0; i < nParam_; ++i) { for (int j = 0; j < nParam_; ++j) { deltaP[i] += alphaRi[i * nParam_ + j] * beta[j]; } }

			for (unsigned short i = 0; i < nParam_; ++i) {
				p_try[i] = p_[i] + deltaP[i];
				p_try[i] = std::min(std::max(pMin_[i], p_try[i]), pMax_[i]);
			}

			for (unsigned short i = 0; i < n_; ++i) {
				deltaY[i] = dat[i] - psf(static_cast<float>(t[i * 2]) - p_try[1], static_cast<float>(t[i * 2 + 1]) - p_try[2],
				                         p_try[0], p_try[3], p_try[4], theta_, p_try[5]);
			}

			// deltaY .* weight
			std::ranges::transform(deltaY, weight, deltaY_w.begin(), std::multiplies<>());
			const float X2_Try = std::transform_reduce(deltaY.begin(), deltaY.end(), deltaY_w.begin(), 0.0f, std::plus<float>(), std::multiplies<float>());

			// lambda * deltaP + beta
			std::ranges::transform(deltaP, beta, beta_l_deltaP.begin(), [lambda](const float d, const float b) { return lambda * d + b; });

			m = std::inner_product(deltaP.begin(), deltaP.end(), beta_l_deltaP.begin(), 0.0f);

			if ((x2 - X2_Try) / (2 * m) > EPSILON_4) { // it IS significantly better
				x2_old = x2;
				std::copy_n(p_try.begin(), nParam_, p_.begin()); // accept p_try

				x2 = lm_matx(t, p_.data(), dat, weight, alpha.data(), beta.data(), hat);
				lambda = std::max(lambda / LAMBDA_DN_FAC, LAMBDA_EPSILON_MIN);
			}
			else {           // it IS NOT better
				x2 = x2_old; // do not accept a_try
				lambda = std::min(lambda * LAMBDA_UP_FAC, LAMBDA_EPSILON_MAX);
			}

			m = 0;
			for (unsigned short i = 0; i < nParam_; ++i) { m = std::max(m, abs(deltaP[i] / p_[i])); }

			if (m < EPSILON_2 && iteration > 2) { stop = true; }
			if (x2 / static_cast<float>(n_) < EPSILON_3 && iteration > 2) { stop = true; }

			m = std::ranges::max(std::views::transform(beta, [](const float b) { return std::abs(b); }));

			if (m < EPSILON_1 && iteration > 2) { stop = true; }
			if (iteration >= iterNum) { stop = true; }
		}

		x = p_[1];
		y = p_[2];
		intensity = p_[0];
		sigmaX = p_[3] / std::numbers::sqrt2_v<float>;
		sigmaY = p_[4] / std::numbers::sqrt2_v<float>;
		offset = p_[5];
		rss = computeRSS(dat, hat);
	}
}
