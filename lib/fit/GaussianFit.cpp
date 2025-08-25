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
	if (roiHalf_ <= Ux && Ux < dim_.x - roiHalf_ && roiHalf_ <= Uy && Uy < dim_.y - roiHalf_) {
		// Création des variables
		double lambda = LAMBDA_0, minI, maxI;
		std::vector<double> alpha(nParams2_), alphaR(nParams2_), alphaRi(nParams2_),
		                    beta(nParams_), beta_l_deltaP(nParams_), p_try(nParams_), deltaP(nParams_),
		                    hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<Pixel> grid(n_);
		std::array<std::vector<double>, 3> p;

		// Initialization of data & t
		preCompute(Ux, Uy, grid, minI, maxI);
		std::ranges::fill(weight, 1);

		initP(p, x, y, z, minI, maxI);

		// Première itération de Levenberg-Marquardt
		double x2     = lm(grid, p[1], weight, alpha, beta, hat);
		double x2_old = x2;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

		bool stop = (m < EPSILON_1);

		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;

			std::ranges::copy(alpha, alphaR.begin());                                                               // Copy du alpha
			for (size_t i = 0; i < nParams_; ++i) { alphaR[i * nParams_ + i] += lambda * alpha[i * nParams_ + i]; } // alphaR = alpha + λ diag(alpha)
			matrix::InverseGaussJordan(alphaR, alphaRi, nParams_);                                                  // Inversion (Gauss-Jordan) et pas LM

			std::ranges::fill(deltaP, 0.0);
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { deltaP[i] += alphaRi[i * nParams_ + j] * beta[j]; }
				p_try[i] = std::clamp(p[1][i] + deltaP[i], p[0][i], p[2][i]);
			}

			// deltaY .* weight
			for (size_t i = 0; i < n_; ++i) { deltaY[i] = grid[i].v - psf(p_try, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y)); }
			std::ranges::transform(deltaY, weight, deltaY_w.begin(), std::multiplies());
			const double X2_Try = std::inner_product(deltaY.begin(), deltaY.end(), deltaY_w.begin(), 0.0);

			// // m = deltaPᵀ * (λ deltaP + beta)
			std::ranges::transform(deltaP, beta, beta_l_deltaP.begin(), [lambda](const double d, const double b) { return lambda * d + b; });
			m = std::inner_product(deltaP.begin(), deltaP.end(), beta_l_deltaP.begin(), 0.0);

			if ((x2 - X2_Try) / (2 * m) > EPSILON_4) { // it IS significantly better
				x2_old = x2;
				p[1]   = p_try;
				x2     = lm(grid, p[1], weight, alpha, beta, hat);
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

		// Mise à jour des différents éléments
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
	const double T1 = std::inner_product(grid.begin(), grid.end(), hat.begin(), 0.0, std::plus(), [](const auto& a, double b) { return a.v * b; }),
	             // data * hat
	             T2 = std::inner_product(grid.begin(), grid.end(), grid.begin(), 0.0, std::plus(),
	                                     [](const auto& a, const auto& b) { return a.v * b.v; }), // data * data
	             T3 = std::inner_product(hat.begin(), hat.end(), hat.begin(), 0.0),               // hat * hat
	             R  = T1 / std::sqrt(T2 * T3);                                                    // Calcul du R

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
double GaussianFit::chi2(std::vector<double>& alpha, std::vector<double>& beta, const std::vector<double>& jacob, const std::vector<double>& weight,
                         const std::vector<double>& delta) const
{
	// tmp = W * J (on multiplie chaque ligne de J par le poids du pixel)
	std::vector<double> tmp(n_ * nParams_);
	for (size_t i = 0; i < n_; ++i) { for (size_t j = 0; j < nParams_; ++j) { tmp[i * nParams_ + j] = weight[i] * jacob[i * nParams_ + j]; } }

	// alpha = J^T * tmp (produit matriciel classique)
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = 0; j < nParams_; ++j) {
			alpha[nParams_ * i + j] = 0;
			for (size_t k = 0; k < n_; ++k) { alpha[nParams_ * i + j] += jacob[k * nParams_ + i] * tmp[k * nParams_ + j]; }
		}
	}

	// beta = J^T * tmp (multiplication matrice-vecteur)
	tmp.resize(n_);
	std::ranges::transform(weight, delta, tmp.begin(), std::multiplies{});
	for (size_t i = 0; i < nParams_; ++i) {
		beta[i] = 0;
		for (size_t j = 0; j < n_; ++j) { beta[i] += jacob[j * nParams_ + i] * tmp[j]; }
	}

	// chi² = ∂ * tmp (somme des carrés pondérée)
	#ifdef __APPLE__
	return std::transform_reduce(delta.begin(), delta.end(), tmp.begin(), 0.0);	 // Pas de support de std::execution sous macOS (Clang/libc++)
	#else
	return std::transform_reduce(std::execution::seq, delta.begin(), delta.end(), tmp.begin(), 0.0);
	#endif
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::lm(const std::vector<Pixel>& grid, const std::vector<double>& p, const std::vector<double>& weight,
                       std::vector<double>& alpha, std::vector<double>& beta, std::vector<double>& hat) const
{
	std::vector<double> delta(n_);           // Différence entre la vraie image (grid[i].v) et la prédiction (hat[i])
	std::vector<double> dydp(n_ * nParams_); //  Jacobienne J (sensibilité de chaque pixel à chaque paramètre)
	// Calcule l'image prédite (hat) et les résidus (delta) pour chaque pixel i de la ROI :
	for (size_t i = 0; i < n_; ++i) {
		hat[i]   = psf(p, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y)); // PSF théorique au pixel (x, y) selon les paramètres p
		delta[i] = grid[i].v - hat[i];                                                     // Résidu = mesure - prédiction
	}

	// Calcule la Jacobienne J (= dydp), c’est-à-dire les dérivées de la PSF par rapport à chacun des paramètres, et ce pour chaque pixel.
	for (size_t i = 0; i < n_; ++i) {
		// Remplit un bloc de nParams_ valeurs dans dydp : [ ∂PSF/∂offset, ∂PSF/∂I0, ∂PSF/∂x0, ∂PSF/∂y0, ... ] au pixel i
		deriv(&dydp[i * nParams_], p, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y));
	}

	// Construit alpha (matrice nParams_ x nParams_), beta (vecteur de taille nParams_), chi² (scalaire).
	// À partir de J (=dydp), des poids W et des résidus r (=delta), ces trois objets sont la "matière première" du pas LM suivant.
	return chi2(alpha, beta, dydp, weight, delta);
}
} // namespace palm_tracer
