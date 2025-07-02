#include "GaussianFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "signal_processing/Matrix.hpp"
#include "signal_processing/PSF.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void GaussianFit::parseParams()
{
	roi_    = static_cast<size_t>(fitParams_[0]);
	sigmaX_ = fitParams_[1];
	sigmaY_ = fitParams_[2];
	theta_  = fitParams_[3];
}

//-------------------------------------------------------------------------------------------------
void GaussianFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point)
{
	double offset = 0, intensity    = 0, rss         = -1,
	       x      = infos.x_, y     = infos.y_, z    = infos.z_,
	       sigmaX = sigmaX_, sigmaY = sigmaY_, theta = theta_;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const auto Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux < width_ - roiHalf_ && roiHalf_ <= Uy && Uy < height_ - roiHalf_) {
		// Création des variables
		double lambda = LAMBDA_0;
		std::vector<double> alpha(nParams2_), alphaR(nParams2_), alphaRi(nParams2_),
		                    beta(nParams_), beta_l_deltaP(nParams_), p_try(nParams_), deltaP(nParams_),
		                    hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<Pixel> grid(n_);

		// Initialization of data & t (A vérifier)
		double minI, maxI;
		preCompute(Ux, Uy, grid, minI, maxI);
		std::ranges::fill(weight, 1);

		std::array<std::vector<double>, 3> p;

		initP(p, x, y, z, minI, maxI);

		double x2     = lm_matx(grid, p[1], weight, alpha, beta, hat);
		double x2_old = x2;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

		bool stop = false;
		if (m < EPSILON_1) { stop = true; }

		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			std::ranges::fill(deltaP, 0.0);

			std::copy_n(alpha.begin(), nParams2_, alphaR.begin()); // Copy du alpha
			for (size_t i = 0; i < nParams_; ++i) { alphaR[i * nParams_ + i] += lambda * alpha[i * nParams_ + i]; }
			InverseGaussJordan(alphaR, alphaRi, nParams_);

			for (size_t i = 0; i < nParams_; ++i) { for (size_t j = 0; j < nParams_; ++j) { deltaP[i] += alphaRi[i * nParams_ + j] * beta[j]; } }

			for (size_t i = 0; i < nParams_; ++i) {
				p_try[i] = p[1][i] + deltaP[i];
				p_try[i] = std::clamp(p_try[i], p[0][i], p[2][i]);
			}

			for (size_t i = 0; i < n_; ++i) { deltaY[i] = grid[i].v - psf(p_try, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y)); }

			// deltaY .* weight
			std::ranges::transform(deltaY, weight, deltaY_w.begin(), std::multiplies());
			const double X2_Try = std::transform_reduce(deltaY.begin(), deltaY.end(), deltaY_w.begin(), 0.0, std::plus(), std::multiplies());

			// lambda * deltaP + beta
			std::ranges::transform(deltaP, beta, beta_l_deltaP.begin(), [lambda](const double d, const double b) { return lambda * d + b; });

			m = std::inner_product(deltaP.begin(), deltaP.end(), beta_l_deltaP.begin(), 0.0);

			if ((x2 - X2_Try) / (2 * m) > EPSILON_4) { // it IS significantly better
				x2_old = x2;
				std::copy_n(p_try.begin(), nParams_, p[1].begin()); // accept p_try

				x2     = lm_matx(grid, p[1], weight, alpha, beta, hat);
				lambda = std::max(lambda / LAMBDA_DN, LAMBDA_MIN);
			}
			else {               // it IS NOT better
				x2     = x2_old; // do not accept a_try
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX);
			}

			m = 0;
			for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::fabs(deltaP[i] / p[1][i])); }

			if (m < EPSILON_2 && iteration > 2) { stop = true; }
			if (x2 / static_cast<double>(n_) < EPSILON_3 && iteration > 2) { stop = true; }

			m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

			if (m < EPSILON_1 && iteration > 2) { stop = true; }
			if (iteration >= N_ITER) { stop = true; }
		}

		offset    = p[1][0];
		intensity = p[1][1];
		x         = p[1][2];
		y         = p[1][3];
		switch (mode_) {
			case fit_mode::none:
			case fit_mode::spline_3d_emccd:
			case fit_mode::spline_3d_scmos:
			case fit_mode::gaussian_2d_fixed: break;
			case fit_mode::gaussian_2d_sigma: {
				sigmaX = p[1][4];
				break;
			}
			case fit_mode::gaussian_3d_fixed: {
				sigmaX = p[1][4] / std::numbers::sqrt2_v<double>;
				sigmaY = p[1][5] / std::numbers::sqrt2_v<double>;
				break;
			}
			case fit_mode::gaussian_3d_theta: {
				sigmaX = p[1][4] / std::numbers::sqrt2_v<double>;
				sigmaY = p[1][5] / std::numbers::sqrt2_v<double>;
				theta  = p[1][6];
				break;
			}
		}
		rss = computeR2(grid, hat);
	}

	switch (mode_) {
		case fit_mode::gaussian_2d_fixed: sigmaX = sigmaY = sigmaX_;
			theta = 0;
			break;
		case fit_mode::gaussian_2d_sigma: sigmaY = sigmaX;
			theta = 0;
			break;
		case fit_mode::none:
		case fit_mode::gaussian_3d_fixed:
		case fit_mode::gaussian_3d_theta:
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: break;
	}
	fill(infos, point, offset, intensity, x, y, z, sigmaX, sigmaY, theta, rss);
}


// ==================== Run Computes ====================
//-------------------------------------------------------------------------------------------------
double GaussianFit::computeR2(std::vector<Pixel>& grid, std::vector<double>& hat) const
{
	// computing "data" mean & "hat" mean
	const double Dat_Mean = std::accumulate(grid.begin(), grid.end(), 0.0,
	                                        [](double sum, const auto& p) { return sum + p.v; }) / static_cast<double>(n_),
	             Hat_Mean = std::accumulate(hat.begin(), hat.end(), 0.0) / static_cast<double>(n_);

	// Normalisation en place
	for (size_t i = 0; i < n_; ++i) {
		grid[i].v -= Dat_Mean;
		hat[i] -= Hat_Mean;
	}

	// Calcul des t1, t2, t3
	const double T1 = std::inner_product(grid.begin(), grid.end(), hat.begin(), 0.0,
	                                     std::plus(), [](const auto& a, double b) { return a.v * b; }), // data * hat
	             T2 = std::inner_product(grid.begin(), grid.end(), grid.begin(), 0.0,
	                                     std::plus(), [](const auto& a, const auto& b) { return a.v * b.v; }), // data * data
	             T3 = std::inner_product(hat.begin(), hat.end(), hat.begin(), 0.0),                            // hat * hat
	             R  = T1 / std::sqrt(T2 * T3);                                                                 // Calcul du R

	return std::isnan(R) ? 0 : R * R;
}

// ==================== PSF ====================
//-------------------------------------------------------------------------------------------------
double GaussianFit::psf(const std::vector<double>& p, const double x, const double y) const
{
	switch (mode_) {
		case fit_mode::gaussian_2d_fixed: return PSF(p[0], p[1], x - p[2], y - p[3], sigmaX_);
		case fit_mode::gaussian_2d_sigma: return PSF(p[0], p[1], x - p[2], y - p[3], p[4]);
		case fit_mode::gaussian_3d_fixed: return PSF(p[0], p[1], x - p[2], y - p[3], p[4], p[5], theta_);
		case fit_mode::gaussian_3d_theta: return PSF(p[0], p[1], x - p[2], y - p[3], p[4], p[5], p[6]);
		case fit_mode::none:
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: return 0;
	}
	return 0;
}


// ==================== Matrix ====================
//-------------------------------------------------------------------------------------------------
void GaussianFit::deriv(double* out, const std::vector<double>& p, const double x, const double y) const
{
	switch (mode_) {
		case fit_mode::none:
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos:
		case fit_mode::gaussian_2d_fixed: Deriv(out, p[1], x - p[2], y - p[3], sigmaX_, false);
			break;
		case fit_mode::gaussian_2d_sigma: Deriv(out, p[1], x - p[2], y - p[3], p[4], true);
			break;
		case fit_mode::gaussian_3d_fixed: Deriv(out, p[1], x - p[2], y - p[3], p[4], p[5], theta_, false);
			break;
		case fit_mode::gaussian_3d_theta: Deriv(out, p[1], x - p[2], y - p[3], p[4], p[5], p[6], true);
			break;
	}
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::chi2(std::vector<double>& alpha, std::vector<double>& beta, const std::vector<double>& dydp, const std::vector<double>& weight,
                         const std::vector<double>& deltaY) const
{
	std::vector<double> tmp(n_ * nParams_);
	for (size_t i = 0; i < n_; ++i) { for (size_t j = 0; j < nParams_; ++j) { tmp[i * nParams_ + j] = weight[i] * dydp[i * nParams_ + j]; } }

	// Computing alpha (A vérifier)
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = 0; j < nParams_; ++j) {
			alpha[nParams_ * i + j] = 0;
			for (size_t k = 0; k < n_; ++k) { alpha[nParams_ * i + j] += dydp[k * nParams_ + i] * tmp[k * nParams_ + j]; }
		}
	}

	// Computing beta
	tmp.resize(n_);
	std::ranges::transform(weight, deltaY, tmp.begin(), std::multiplies{});

	for (size_t i = 0; i < nParams_; ++i) {
		beta[i] = 0;
		for (size_t j = 0; j < n_; ++j) { beta[i] += dydp[j * nParams_ + i] * tmp[j]; }
	}

	// Computing Chi_sq
	#ifdef __APPLE__
		// Pas de support de std::execution sous macOS (Clang/libc++)
		return std::transform_reduce(deltaY.begin(), deltaY.end(), tmp.begin(), 0.0);
	#else
	return std::transform_reduce(std::execution::seq, deltaY.begin(), deltaY.end(), tmp.begin(), 0.0);
	#endif
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::lm_matx(const std::vector<Pixel>& grid, const std::vector<double>& p, const std::vector<double>& weight,
                            std::vector<double>& alpha, std::vector<double>& beta, std::vector<double>& hat) const
{
	std::vector<double> delta(n_);
	for (size_t i = 0; i < n_; ++i) {
		hat[i]   = psf(p, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y));
		delta[i] = grid[i].v - hat[i];
	}

	// Computing Jacobian matrix (dydp)
	std::vector<double> dydp(n_ * nParams_);
	for (size_t i = 0; i < n_; ++i) { deriv(&dydp[i * nParams_], p, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y)); }

	return chi2(alpha, beta, dydp, weight, delta);
}
} // namespace palm_tracer
