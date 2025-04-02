#include "GaussianFit2D.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <ctime>

float GaussianFit2D::lm_matx(const std::vector<unsigned int>& t, const float* p, const std::vector<float>& dat, const std::vector<float>& weight, float* alpha, float* beta, std::vector<float>& hat)
{
	std::vector<float> delta_y(n_);
	for (size_t i = 0; i < n_; ++i) {
		hat[i] = psf(static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], sigma_, p[3]);
		delta_y[i] = dat[i] - hat[i];
	}

	// Computing Jacobian matrix (dydp)
	std::vector<float> dydp(n_ * nParam_);
	for (size_t i = 0; i < n_; ++i) {
		deriv(&dydp[i * nParam_], static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], sigma_, false);
	}

	return chi2(alpha, beta, dydp, weight, delta_y);
}

void GaussianFit2D::initP(const float x, const float y)
{
	pMin_ = {0, x - 1, y - 1, 0};
	p_ = {1000, x, y, 1000};
	pMax_ = {65535, x + 1, y + 1, 65535};
}


void GaussianFit2D::Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset,
                        float& RSS, int iterNum, unsigned short size)
{
	if (size % 2 == 0) { size--; } // Not present in GPU original file why ?

	n_ = size * size;
	const unsigned short Half = size / 2;
	sigma_ = sigmaX;

	float X2;
	if (x0 >= Half && x0 <= sizeX - Half - 1 && y0 >= Half && y0 <= sizeY - Half - 1) {
		initP(x0, y0);

		std::vector<float> dat(n_), hat(n_), weight(n_);
		std::vector<unsigned int> t(2 * n_);

		// Initialization of dat & t (A vérifier)
		float minI, maxI; // Not present in GPU original file why ?
		preComputeFit(minI, maxI, image, sizeX, static_cast<int>(x0), static_cast<int>(y0), dat, weight, t, size);
		p_[0] = maxI - minI; // Not present in GPU original file why ?
		p_[3] = minI;        // Not present in GPU original file why ?

		float lambda;
		std::array<float, N_PARAM2> alpha;
		std::array<float, N_PARAM> beta;

		float X2_old;
		X2 = lm_matx(t, p_.data(), dat, weight, alpha.data(), beta.data(), hat);

		float max_ = 0;
		for (int i = 0; i < nParam_; ++i) { max_ = std::max(max_, abs(beta[i])); }

		unsigned int stop = 0;
		if (max_ < EPSILON_1) { stop = 1; }

		lambda = LAMBDA_0;
		X2_old = X2;

		unsigned short iteration = 0;

		std::array<float, N_PARAM2> inter;
		std::array<float, N_PARAM2> inter1;
		std::vector<float> inter2(n_);
		std::array<float, N_PARAM> inter3;

		std::array<float, N_PARAM> p_try;
		std::vector<float> delta_y(n_);
		float X2_try;
		float m;
		float rho;
		std::array<float, N_PARAM> delta_p;  // Déclaré une seule fois avant la boucle

		while (stop == 0 && iteration < iterNum) // Main Loop
		{
			iteration++;
			delta_p.fill(0.0f);

			// Copy alpha to inter
			for (int i = 0; i < nParam2_; ++i) { inter[i] = alpha[i]; }
			for (int i = 0; i < nParam_; ++i) { inter[i * nParam_ + i] += lambda * alpha[i * nParam_ + i]; }

			inverse(inter.data(), inter1.data());

			for (int i = 0; i < nParam_; ++i) {
				for (int cursor = 0; cursor < nParam_; cursor++) { delta_p[i] += inter1[i * nParam_ + cursor] * beta[cursor]; }
			}

			for (int i = 0; i < nParam_; ++i) {
				p_try[i] = p_[i] + delta_p[i];
				p_try[i] = std::min(std::max(pMin_[i], p_try[i]), pMax_[i]);
			}

			for (int i = 0; i < n_; ++i) {
				delta_y[i] = dat[i] - psf((float)t[i * 2] - p_try[1], (float)t[i * 2 + 1] - p_try[2], p_try[0], sigma_, p_try[3]);
			}

			for (int i = 0; i < n_; ++i) { inter2[i] = delta_y[i] * weight[i]; } // delta_y .* weight

			X2_try = 0;
			for (int i = 0; i < n_; ++i) { X2_try += delta_y[i] * inter2[i]; }
			for (int i = 0; i < nParam_; ++i) { inter3[i] = lambda * delta_p[i] + beta[i]; } // lambda * delta_p + beta

			m = 0;
			for (int i = 0; i < nParam_; ++i) { m += delta_p[i] * inter3[i]; }
			m *= 2;

			rho = (X2 - X2_try) / m;

			if (rho > EPSILON_4) // it IS significantly better
			{
				X2_old = X2;
				for (int i = 0; i < nParam_; ++i) { p_[i] = p_try[i]; } // accept p_try
				X2 = lm_matx(t, p_.data(), dat, weight, alpha.data(), beta.data(), hat);
				lambda = std::max(lambda / LAMBDA_DN_FAC, 1.e-7f);
			}
			else // it IS NOT better
			{
				X2 = X2_old; // do not accept a_try
				lambda = std::min(lambda * LAMBDA_UP_FAC, 1.e7f);
			}

			max_ = 0;
			for (int i = 0; i < nParam_; ++i) { max_ = std::max(max_, abs(delta_p[i] / p_[i])); }
			if (max_ < EPSILON_2 && iteration > 2) { stop = 1; }
			if (X2 / n_ < EPSILON_3 && iteration > 2) { stop = 1; }

			max_ = 0;
			for (int i = 0; i < nParam_; ++i) { max_ = std::max(max_, abs(beta[i])); }
			if (max_ < EPSILON_1 && iteration > 2) { stop = 1; }
			if (iteration == iterNum) { stop = 1; }
		}

		x0 = p_[1];
		y0 = p_[2];
		I0 = p_[0];
		offset = p_[3];
		RSS = computeRSS(dat, hat);
	}
}
