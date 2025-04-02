#include "GaussianFit3D.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numbers>

float GaussianFit3D::lm_matx(std::vector<unsigned int>& t, float* p, std::vector<float>& y_dat, std::vector<float>& weight, float* alpha, float* beta, std::vector<float>& y_hat)
{
	std::vector<float> delta_y(n_);
	for (size_t i = 0; i < n_; ++i) {
		y_hat[i] = psf(static_cast<float>(t[i * 2]) - p[1], static_cast<float>(t[i * 2 + 1]) - p[2], p[0], p[3], p[4], theta_, p[5]);
		delta_y[i] = y_dat[i] - y_hat[i];
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

void GaussianFit3D::Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& sigmaY,
                        float& angle, float& offset, float& RSS, int iterNum, unsigned short size)
{
	if (size % 2 == 0) { size--; } // Not present in GPU original file why ?

	n_ = size * size;
	const unsigned short Half = size / 2;
	theta_ = angle;

	float X2;
	if (x0 >= (int)Half && x0 <= sizeX - (int)Half - 1 && y0 >= (int)Half && y0 <= sizeY - (int)Half - 1) {
		// for real data
		initP(x0, y0, sigmaX, sigmaY);

		std::vector<float> y_dat(n_), y_hat(n_), weight(n_);
		std::vector<unsigned int> t(2 * n_);

		// Initialization of y_dat & t (A vérifier)
		float minI, maxI; // Not present in GPU original file why ?
		preComputeFit(minI, maxI, image, sizeX, static_cast<int>(x0), static_cast<int>(y0), y_dat, weight, t, size);
		p_[0] = maxI - minI; // Not present in GPU original file why ?
		p_[5] = minI;        // Not present in GPU original file why ?

		float lambda;
		float alpha[N_PARAM2];
		float beta[N_PARAM];

		float X2_old;
		X2 = lm_matx(t, p_.data(), y_dat, weight, alpha, beta, y_hat);

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

		while (stop == 0 && iteration <= iterNum) // Main Loop
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
				delta_y[i] = y_dat[i] - psf((float)t[i * 2] - p_try[1], (float)t[i * 2 + 1] - p_try[2], p_try[0], p_try[3], p_try[4], theta_, p_try[5]);
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

				X2 = lm_matx(t, p_.data(), y_dat, weight, alpha, beta, y_hat);

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

		// computing "y_dat" mean & "y_hat" mean
		float y_dat_mean = 0.;
		float y_hat_mean = 0.;
		for (int i = 0; i < n_; ++i) {
			y_dat_mean += y_dat[i];
			y_hat_mean += y_hat[i];
		}
		y_dat_mean /= n_;
		y_hat_mean /= n_;

		// computing 't1' & 't2' & 't3'
		float t1 = 0.;
		float t2 = 0.;
		float t3 = 0.;
		float diff_hat, diff_dat;

		for (int i = 0; i < n_; ++i) {
			diff_dat = y_dat[i] - y_dat_mean;
			diff_hat = y_hat[i] - y_hat_mean;
			t1 += diff_dat * diff_hat;
			t2 += diff_dat * diff_dat;
			t3 += diff_hat * diff_hat;
		}
		t2 = sqrt(t2);
		t3 = sqrt(t3);

		// computing 'R'
		const float R = t1 / (t2 * t3);

		x0 = p_[1];
		y0 = p_[2];
		I0 = p_[0];
		sigmaX = p_[3] / sqrt(2.);
		sigmaY = p_[4] / sqrt(2.);
		offset = p_[5];

		if (_isnan(R)) { RSS = static_cast<float>(0); }
		else { RSS = R * R; }
	}
}
