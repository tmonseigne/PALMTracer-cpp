#include "GaussianFit2D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <ranges>

double GaussianFit2D::lm_matx(const std::vector<size_t>& t, const double* p, const std::vector<double>& dat, const std::vector<double>& weight,
                              double* alpha, double* beta, std::vector<double>& hat)
{
	std::vector<double> delta(n_);
	for (size_t i = 0; i < n_; ++i) {
		hat[i] = psf(static_cast<double>(t[i * 2]) - p[1], static_cast<double>(t[i * 2 + 1]) - p[2], p[0], sigma_, p[3]);
		delta[i] = dat[i] - hat[i];
	}

	// Computing Jacobian matrix (dydp)
	std::vector<double> dydp(n_ * nParam_);
	for (size_t i = 0; i < n_; ++i) {
		deriv(&dydp[i * nParam_], static_cast<double>(t[i * 2]) - p[1], static_cast<double>(t[i * 2 + 1]) - p[2], p[0], sigma_, false);
	}

	return chi2(alpha, beta, dydp, weight, delta);
}

void GaussianFit2D::initP(const double x, const double y)
{
	pMin_ = {0, x - 1, y - 1, 0};
	p_ = {1000, x, y, 1000};
	pMax_ = {65535, x + 1, y + 1, 65535};
}


void GaussianFit2D::Run(const unsigned short* image, const size_t sizeX, const size_t sizeY,
                        double& x, double& y, double& intensity, const double& sigma, double& offset,
                        double& rss, const size_t iterNum, size_t roi)
{
	if (roi % 2 == 0) { roi--; } // Not present in GPU original file why ?
	n_ = roi * roi;
	sigma_ = sigma;

	const size_t Half = roi / 2, Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (Half <= Ux && Ux <= sizeX - Half - 1 && Half <= Uy && Uy <= sizeY - Half - 1) {
		initP(x, y);

		// Création des variables
		double lambda = LAMBDA_0;

		std::array<double, N_PARAM2> alpha, alphaR, alphaRi;
		std::array<double, N_PARAM> beta, beta_l_deltaP, p_try, deltaP;
		std::vector<double> dat(n_), hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<size_t> t(2 * n_);

		// Initialization of dat & t (A vérifier)
		double minI, maxI; // Not present in GPU original file why ?
		preComputeFit(minI, maxI, image, sizeX, static_cast<int>(x), static_cast<int>(y), dat, weight, t, roi);
		p_[0] = maxI - minI; // Not present in GPU original file why ?
		p_[3] = minI;        // Not present in GPU original file why ?

		double x2 = lm_matx(t, p_.data(), dat, weight, alpha.data(), beta.data(), hat);
		double x2_old = x2;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::abs(b); }));

		bool stop = false;
		if (m < EPSILON_1) { stop = true; }

		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			deltaP.fill(0.0);

			std::copy_n(alpha.begin(), nParam2_, alphaR.begin()); // Copy du alpha
			for (size_t i = 0; i < nParam_; ++i) { alphaR[i * nParam_ + i] += lambda * alpha[i * nParam_ + i]; }
			inverse(alphaR.data(), alphaRi.data());

			for (size_t i = 0; i < nParam_; ++i) { for (size_t j = 0; j < nParam_; ++j) { deltaP[i] += alphaRi[i * nParam_ + j] * beta[j]; } }

			for (size_t i = 0; i < nParam_; ++i) {
				p_try[i] = p_[i] + deltaP[i];
				p_try[i] = std::min(std::max(pMin_[i], p_try[i]), pMax_[i]);
			}

			for (size_t i = 0; i < n_; ++i) {
				deltaY[i] = dat[i] - psf(static_cast<double>(t[i * 2]) - p_try[1], static_cast<double>(t[i * 2 + 1]) - p_try[2], p_try[0], sigma_, p_try[3]);
			}

			// deltaY .* weight
			std::ranges::transform(deltaY, weight, deltaY_w.begin(), std::multiplies<>());
			const double X2_Try = std::transform_reduce(deltaY.begin(), deltaY.end(), deltaY_w.begin(), 0.0, std::plus<double>(), std::multiplies<double>());

			// lambda * deltaP + beta
			std::ranges::transform(deltaP, beta, beta_l_deltaP.begin(), [lambda](const double d, const double b) { return lambda * d + b; });

			m = std::inner_product(deltaP.begin(), deltaP.end(), beta_l_deltaP.begin(), 0.0);

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
			for (size_t i = 0; i < nParam_; ++i) { m = std::max(m, abs(deltaP[i] / p_[i])); }

			if (m < EPSILON_2 && iteration > 2) { stop = true; }
			if (x2 / static_cast<double>(n_) < EPSILON_3 && iteration > 2) { stop = true; }

			m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::abs(b); }));

			if (m < EPSILON_1 && iteration > 2) { stop = true; }
			if (iteration >= iterNum) { stop = true; }
		}

		x = p_[1];
		y = p_[2];
		intensity = p_[0];
		offset = p_[3];
		rss = computeRSS(dat, hat);
	}
}
