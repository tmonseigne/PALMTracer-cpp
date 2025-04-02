#include "GaussianFit2D.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <ctime>

float GaussianFit2D::lm_matx(unsigned int* t, float* p, float* y_dat, float* weight_sq, float* alpha, float* beta, float* y_hat)
{
	float* delta_y = (float*)malloc(n_ * sizeof(float));
	for (int elt = 0; elt < n_; elt++) {
		y_hat[elt] = psf((float)t[elt * 2] - p[1], (float)t[elt * 2 + 1] - p[2], p[0], sigma_, p[3]);
		delta_y[elt] = y_dat[elt] - y_hat[elt];
	}

	// Computing Jacobian matrix (dydp)
	float* dydp = (float*)malloc(n_ * nParam_ * sizeof(float));
	for (size_t elt = 0; elt < n_; elt++) { deriv(&dydp[elt * nParam_], (float)t[elt * 2] - p[1], (float)t[elt * 2 + 1] - p[2], p[0], sigma_, false); }

	const float Res = chi2(alpha, beta, dydp, weight_sq, delta_y);

	// cleanup memory
	free(delta_y);
	free(dydp);

	return Res;
}

void GaussianFit2D::Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset,
                        float& RSS, int iterNum, unsigned short size)
{
	if (size % 2 == 0) {
		size--; // Not present in GPU original file why ?
	}

	dimFit_ = size;
	n_ = dimFit_ * dimFit_;
	ref_ = dimFit_ / 2;
	sigma_ = sigmaX;

	float X2;

	if (x0 >= (int)ref_ && x0 <= sizeX - (int)ref_ - 1 && y0 >= (int)ref_ && y0 <= sizeY - (int)ref_ - 1) {
		// for real data
		float p_min[4] = {0., x0 - 1, y0 - 1, 0.};
		float p[4] = {1000., x0, y0, 1000.};
		float p_max[4] = {65535., x0 + 1, y0 + 1, 65535.};

		float* y_dat = (float*)malloc(n_ * sizeof(float));
		float* weight = (float*)malloc(n_ * sizeof(float));
		float* weight_sq = (float*)malloc(n_ * sizeof(float));
		unsigned int* t = (unsigned int*)malloc(2 * n_ * sizeof(unsigned int));

		// Initialization of y_dat & t (A vérifier)
		float max_intensity = -70000; // Not present in GPU original file why ?
		float min_intensity = +70000; // Not present in GPU original file why ?
		for (int jj = -dimFit_ / 2; jj <= dimFit_ / 2; jj++) {
			for (int ii = -dimFit_ / 2; ii <= dimFit_ / 2; ii++) {
				y_dat[ii + dimFit_ / 2 + (jj + dimFit_ / 2) * dimFit_] = (float)image[(int)x0 + ii + ((int)y0 + jj) * sizeX];
				t[(ii + dimFit_ / 2 + (jj + dimFit_ / 2) * dimFit_) * 2] = (int)x0 + ii;
				t[(ii + dimFit_ / 2 + (jj + dimFit_ / 2) * dimFit_) * 2 + 1] = (int)y0 + jj;

				max_intensity = std::max(max_intensity, y_dat[ii + dimFit_ / 2 + (jj + dimFit_ / 2) * dimFit_]); // Not present in GPU original file why ?
				min_intensity = std::min(min_intensity, y_dat[ii + dimFit_ / 2 + (jj + dimFit_ / 2) * dimFit_]); // Not present in GPU original file why ?
			}
		}

		p[0] = max_intensity - min_intensity; // Not present in GPU original file why ?
		p[3] = min_intensity;                 // Not present in GPU original file why ?

		for (int ii = 0; ii < n_; ii++) {
			weight[ii] = 1/*sqrt((n_-nParam_+1) / (y_dat[ii]*y_dat[ii]))*/;
			weight_sq[ii] = pow(weight[ii], 2);
		}

		constexpr float epsilon_1 = 1e-6;
		constexpr float epsilon_2 = 1e-6;
		constexpr float epsilon_3 = 1e-5;
		constexpr float epsilon_4 = 1e-2;
		constexpr float lambda_0 = 1e-2;
		constexpr float lambda_UP_fac = 11;
		constexpr float lambda_DN_fac = 9;

		float lambda;

		float alpha[16];
		float beta[4];
		//float y_hat[n_];
		float* y_hat = (float*)malloc(n_ * sizeof(float));

		//float X2;
		float X2_old;
		X2 = lm_matx(t, p, y_dat, weight_sq, alpha, beta, y_hat);

		float max_ = 0;
		for (int ii = 0; ii < nParam_; ii++) { max_ = std::max(max_, abs(beta[ii])); }

		unsigned int stop = 0;
		if (max_ < epsilon_1) { stop = 1; }

		lambda = lambda_0;
		X2_old = X2;

		unsigned short iteration = 0;

		float inter[16];
		float inter1[16];
		float* inter2 = (float*)malloc(n_ * sizeof(float));
		float inter3[4];

		float p_try[4];
		float* delta_y = (float*)malloc(n_ * sizeof(float));
		float X2_try;
		float m;
		float rho;

		while (stop == 0 && iteration < iterNum) // Main Loop
		{
			iteration++;
			float delta_p[4] = {0, 0, 0, 0};

			// Copy alpha to inter
			for (int ii = 0; ii < nParam_ * nParam_; ii++) { inter[ii] = alpha[ii]; }

			inter[0 * nParam_ + 0] += lambda * alpha[0 * nParam_ + 0];
			inter[1 * nParam_ + 1] += lambda * alpha[1 * nParam_ + 1];
			inter[2 * nParam_ + 2] += lambda * alpha[2 * nParam_ + 2];
			inter[3 * nParam_ + 3] += lambda * alpha[3 * nParam_ + 3];

			inverse(inter, inter1);

			for (int ii = 0; ii < nParam_; ii++) {
				for (int cursor = 0; cursor < nParam_; cursor++) { delta_p[ii] += inter1[ii * nParam_ + cursor] * beta[cursor]; }
			}

			for (int ii = 0; ii < nParam_; ii++) {
				p_try[ii] = p[ii] + delta_p[ii];
				p_try[ii] = std::min(std::max(p_min[ii], p_try[ii]), p_max[ii]);
			}

			for (int ii = 0; ii < n_; ii++) {
				delta_y[ii] = y_dat[ii] - psf((float)t[ii * 2] - p_try[1], (float)t[ii * 2 + 1] - p_try[2], p_try[0], sigma_, p_try[3]);
			}

			for (int ii = 0; ii < n_; ii++) { inter2[ii] = delta_y[ii] * weight_sq[ii]; } // delta_y .* weight_sq

			X2_try = 0;
			for (int ii = 0; ii < n_; ii++) { X2_try += delta_y[ii] * inter2[ii]; }

			for (int ii = 0; ii < nParam_; ii++) { inter3[ii] = lambda * delta_p[ii] + beta[ii]; } // lambda * delta_p + beta

			m = 0;
			for (int ii = 0; ii < nParam_; ii++) { m += delta_p[ii] * inter3[ii]; }
			m *= 2;

			rho = (X2 - X2_try) / m;

			if (rho > epsilon_4) // it IS significantly better
			{
				X2_old = X2;
				for (int ii = 0; ii < nParam_; ii++) { p[ii] = p_try[ii]; } // accept p_try

				X2 = lm_matx(t, p, y_dat, weight_sq, alpha, beta, y_hat);

				lambda = std::max(lambda / lambda_DN_fac, 1.e-7f);
			}
			else // it IS NOT better
			{
				X2 = X2_old; // do not accept a_try
				lambda = std::min(lambda * lambda_UP_fac, 1.e7f);
			}

			max_ = 0;
			for (int ii = 0; ii < nParam_; ii++) { max_ = std::max(max_, abs(delta_p[ii] / p[ii])); }

			if (max_ < epsilon_2 && iteration > 2) { stop = 1; }

			if (X2 / n_ < epsilon_3 && iteration > 2) { stop = 1; }

			max_ = 0;
			for (int ii = 0; ii < nParam_; ii++) { max_ = std::max(max_, abs(beta[ii])); }

			if (max_ < epsilon_1 && iteration > 2) { stop = 1; }

			if (iteration == iterNum) { stop = 1; }
		}

		// Computing statistics

		// computing "y_dat" mean & "y_hat" mean
		float y_dat_mean = 0.;
		float y_hat_mean = 0.;
		for (int ii = 0; ii < n_; ii++) {
			y_dat_mean += y_dat[ii];
			y_hat_mean += y_hat[ii];
		}
		y_dat_mean /= float(n_);
		y_hat_mean /= float(n_);

		// computing 't1' & 't2' & 't3'
		float t1 = 0.;
		float t2 = 0.;
		float t3 = 0.;
		float diff_hat, diff_dat;

		for (int ii = 0; ii < n_; ii++) {
			diff_dat = y_dat[ii] - y_dat_mean;
			diff_hat = y_hat[ii] - y_hat_mean;
			t1 += diff_dat * diff_hat;
			t2 += diff_dat * diff_dat;
			t3 += diff_hat * diff_hat;
		}
		t2 = sqrt(t2);
		t3 = sqrt(t3);

		// computing 'R_sq'
		float R_sq = pow(t1 / (t2 * t3), 2);

		x0 = p[1];
		y0 = p[2];
		I0 = p[0];
		offset = p[3];
		if (_isnan(R_sq)) { RSS = static_cast<float>(0); }
		else { RSS = R_sq; }

		// cleanup memory
		free(y_dat);
		free(weight);
		free(weight_sq);
		free(t);
		free(y_hat);
		free(inter2);
		free(delta_y);
	}
}
