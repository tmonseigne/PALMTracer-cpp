#include "GaussianFit3D_Theta.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>

float GaussianFit3D_Theta::lm_matx(const std::vector<unsigned int>& t, const float* p, const std::vector<float>& dat, const std::vector<float>& weight,
                                   float* alpha, float* beta, std::vector<float>& hat)
{
	std::vector<float> delta(n_);
	for (size_t i = 0; i < n_; ++i) {
		hat[i] = psf(static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], p[3], p[4], p[5], p[6]);
		delta[i] = dat[i] - hat[i];
	}

	// Computing Jacobian matrix (dydp)
	std::vector<float> dydp(static_cast<size_t>(n_ * nParam_));
	for (size_t i = 0; i < n_; ++i) {
		deriv(&dydp[i * nParam_], static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], p[3], p[4], p[5], true);
	}

	return chi2(alpha, beta, dydp, weight, delta);
}

void GaussianFit3D_Theta::initP(const float x, const float y, const float sigmaX, const float sigmaY, const float theta)
{
	pMin_ = {0, x - 1, y - 1, 0.1f, 0.1f, -150, 0};
	p_ = {2000, x, y, std::numbers::sqrt2_v<float> * sigmaX, std::numbers::sqrt2_v<float> * sigmaY, theta, 700};
	pMax_ = {65535, x + 1, y + 1, 20, 20, 150, 65535};
}

void GaussianFit3D_Theta::Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& i0, float& sigmaX, float& sigmaY,
                              float& angle, float& offset, float& rss, int iterNum, unsigned short size)
{
	if (size % 2 == 0) { size--; } // Not present in GPU original file why ?

	n_ = size * size;
	const unsigned short Half = size / 2;

	if (Half <= x0 && x0 <= sizeX - Half - 1 && Half <= y0 && y0 <= sizeY - Half - 1) {
		initP(x0, y0, sigmaX, sigmaY, angle);

		// Création des variables
		float lambda = LAMBDA_0;

		std::array<float, N_PARAM2> alpha, alphaR, alphaRi;
		std::array<float, N_PARAM> beta, beta_l_deltaP, p_try, deltaP;
		std::vector<float> dat(n_), hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<unsigned int> t(2 * n_);

		// Initialization of dat & t (A vérifier)
		float minI, maxI; // Not present in GPU original file why ?
		preComputeFit(minI, maxI, image, sizeX, static_cast<int>(x0), static_cast<int>(y0), dat, weight, t, size);
		p_[0] = maxI - minI; // Not present in GPU original file why ?
		p_[6] = minI;        // Not present in GPU original file why ?

		float x2 = lm_matx(t, p_.data(), dat, weight, alpha.data(), beta.data(), hat);
		float x2_old = x2;

		// Trouver la valeur absolue maximale
		float m = std::ranges::max(std::views::transform(beta, [](const float x) { return std::abs(x); }));

		unsigned int stop = 0;
		if (m < EPSILON_1) { stop = 1; }


		unsigned short iteration = 0;
		while (stop == 0 && iteration <= iterNum) // Main Loop
		{
			iteration++;
			deltaP.fill(0.0f);

			std::copy_n(alpha.begin(), nParam2_, alphaR.begin()); // Copy du alpha
			for (unsigned short i = 0; i < nParam_; ++i) { alphaR[i * nParam_ + i] += lambda * alpha[i * nParam_ + i]; }

			inverse(alphaR.data(), alphaRi.data());

			for (unsigned short i = 0; i < nParam_; ++i) {
				for (int j = 0; j < nParam_; ++j) { deltaP[i] += alphaRi[i * nParam_ + j] * beta[j]; }
			}

			for (unsigned short i = 0; i < nParam_; ++i) {
				p_try[i] = p_[i] + deltaP[i];
				p_try[i] = std::min(std::max(pMin_[i], p_try[i]), pMax_[i]);
			}

			for (unsigned short i = 0; i < n_; ++i) {
				deltaY[i] = dat[i] - psf(static_cast<float>(t[i * 2]) - p_try[1], static_cast<float>(t[i * 2 + 1]) - p_try[2],
				                         p_try[0], p_try[3], p_try[4], p_try[5], p_try[6]);
			}

			// deltaY .* weight
			std::transform(deltaY.begin(), deltaY.end(), weight.begin(), deltaY_w.begin(), std::multiplies<float>());
			const float X2_Try = std::transform_reduce(deltaY.begin(), deltaY.end(), deltaY_w.begin(), 0.0f, std::plus<float>(), std::multiplies<float>());

			// lambda * deltaP + beta
			std::transform(deltaP.begin(), deltaP.end(), beta.begin(), beta_l_deltaP.begin(), [lambda](float delta, float b) { return lambda * delta + b; });

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

			if (m < EPSILON_2 && iteration > 2) { stop = 1; }
			if (x2 / static_cast<float>(n_) < EPSILON_3 && iteration > 2) { stop = 1; }

			m = std::ranges::max(std::views::transform(beta, [](const float x) { return std::abs(x); }));

			if (m < EPSILON_1 && iteration > 2) { stop = 1; }
			if (iteration == iterNum) { stop = 1; }
		}

		x0 = p_[1];
		y0 = p_[2];
		i0 = p_[0];
		sigmaX = p_[3] / std::numbers::sqrt2_v<float>;
		sigmaY = p_[4] / std::numbers::sqrt2_v<float>;
		angle = p_[5];
		offset = p_[6]/*iteration*/;
		rss = computeRSS(dat, hat);
	}
}
