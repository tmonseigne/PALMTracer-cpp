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
void GaussianFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const
{
	// Valeurs par défaut
	double offset = 0, intensity    = 0, rss         = -1,
	       x      = infos.x_, y     = infos.y_, z    = infos.z_,
	       sigmaX = sigmaX_, sigmaY = sigmaY_, theta = theta_;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const auto Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux < dim_.x - roiHalf_ && roiHalf_ <= Uy && Uy < dim_.y - roiHalf_) {
		// Création des variables
		double lambda = LAMBDA_0, minI = 0, maxI = 0;
		std::vector<double> alpha(nParams2_), alphaR(nParams2_), alphaRi(nParams2_),
		                    beta(nParams_), grad(nParams_), deltaP(nParams_), lBetaDeltaP(nParams_), p_try(nParams_),
		                    hat(n_), weight(n_, 1), delta(n_), deltaW(n_);
		std::vector<Pixel> grid(n_);
		std::array<std::vector<double>, 3> p;

		// Initialisation des données (Récupération des min/max, de la ROI et initialisation du P et de ses bornes)
		preCompute(Ux, Uy, grid, minI, maxI);
		initP(p, x, y, z, minI, maxI);

		// Première itération de Levenberg-Marquardt
		double chi2 = lm(grid, p[1], weight, hat, alpha, beta, delta, deltaW, grad);
		double chi2_old = chi2;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

		const bool Stop = (m < EPSILON_1);

		for (size_t iteration = (Stop ? N_ITER : 0); iteration < N_ITER; ++iteration) { // Main Loop
			alphaR = alpha;
			for (size_t i = 0; i < nParams_; ++i) { alphaR[i * nParams_ + i] += lambda * alpha[i * nParams_ + i]; } // alphaR = alpha + λ diag(alpha)
			matrix::InverseGaussJordan(alphaR, alphaRi, nParams_);                                                  // Inversion (Gauss-Jordan) et pas LM

			std::ranges::fill(deltaP, 0.0);
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { deltaP[i] += alphaRi[i * nParams_ + j] * beta[j]; } // δP
				p_try[i] = std::clamp(p[1][i] + deltaP[i], p[0][i], p[2][i]);                               // p_try = clamp(p + δP)
			}

			// Gain prédit (pas d’éval d’image ici)
			// m = deltaPᵀ * (λ deltaP + beta)
			for (size_t i = 0; i < nParams_; ++i) { lBetaDeltaP[i] = lambda * deltaP[i] + beta[i]; }
			m = std::inner_product(deltaP.begin(), deltaP.end(), lBetaDeltaP.begin(), 0.0);

			const double Gain = gain(deltaP, alpha, beta);
			if (Gain <= 0.0) {
				// clairement pas bon → on n’évalue pas, on augmente λ et on retente
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de lambda
				// critères d’arrêt “petit pas”
				m = 0.0;
				for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::fabs(deltaP[i]) / std::max(std::fabs(p[1][i]), 1e-12)); }
				if (m < EPSILON_2 && iteration >= 2) { break; } // On est trop faible est déjà 2 iterations, on arrête
				continue;                                       // On repart au début de la boucle
			}

			// On calcule les nouvelles valeurs : évalue image + lm pour p_try
			double chi2_try = 0.0;
			for (size_t i = 0; i < n_; ++i) {
				const double Dx = grid[i].v - psf(p_try, static_cast<double>(grid[i].x), static_cast<double>(grid[i].y));
				chi2_try += weight[i] * Dx * Dx;
			}
			const double Rho = (chi2 - chi2_try) / (2.0 * m);

			if (Rho > EPSILON_4) { // On accepte : on remplace l’état courant par _try
				// On accepte : on remplace l’état courant par _try
				chi2_old = chi2;
				p[1] = p_try;
				chi2 = lm(grid, p[1], weight, hat, alpha, beta, delta, deltaW, grad);
				lambda = std::max(lambda / LAMBDA_DN, LAMBDA_MIN);
			}
			else {                                                 // On refuse : on garde l’état courant
				chi2 = chi2_old;
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de lambda
			}

			// Critères d'arrêts
			m = 0.0;
			for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::abs(deltaP[i]) / std::max(std::abs(p[1][i]), 1e-12)); }
			if (iteration >= 2) { // Si on est déjà à 2 iterations, on arrête
				if (m < EPSILON_2) { break; }
				if (chi2 / static_cast<double>(n_) < EPSILON_3) { break; }
				if (std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); })) < EPSILON_1) { break; }
			}
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
		rss = computeR2(grid, hat); // Calcul du RSS sur XY
	}

	// Post-normalisation selon mode
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
double GaussianFit::lm(const std::vector<Pixel>& grid, const std::vector<double>& p, const std::vector<double>& weight, std::vector<double>& hat,
                       std::vector<double>& alpha, std::vector<double>& beta, std::vector<double>& delta, std::vector<double>& deltaW,
                       std::vector<double>& grad) const
{
	std::ranges::fill(alpha, 0.0);
	std::ranges::fill(beta, 0.0);
	double chi2 = 0.0;

	for (size_t i = 0; i < n_; ++i) {
		const double X = static_cast<double>(grid[i].x), Y = static_cast<double>(grid[i].y);
		hat[i]         = psf(p, X, Y);         // PSF théorique au pixel (x, y) selon les paramètres p
		delta[i]       = grid[i].v - hat[i];   // Résidu (∂) = mesure - prédiction
		deltaW[i]      = weight[i] * delta[i]; // Résidu pondéré (W * ∂)
		chi2 += delta[i] * deltaW[i];          // mise à jour du Chi² (somme des carrés pondérée W * ∂ * ∂)
		deriv(grad.data(), p, X, Y);           // dérivées au pixel

		// Accumulations alpha += w * grad * gradᵀ  et  beta += w * grad * r
		for (size_t a = 0; a < nParams_; ++a) {
			beta[a] += grad[a] * deltaW[i];
			// rang-1 symétrique
			const size_t Row = a * nParams_;
			for (size_t b = 0; b < nParams_; ++b) { alpha[Row + b] += weight[i] * grad[a] * grad[b]; }
		}
	}
	return chi2;
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::gain(const std::vector<double>& deltaP, const std::vector<double>& alpha, const std::vector<double>& beta) const
{
	// t = alpha * deltaP
	double quad = 0.0, lin = 0.0;
	for (size_t i = 0; i < nParams_; ++i) {
		double ti        = 0.0;
		const size_t Row = i * nParams_;
		for (size_t j = 0; j < nParams_; ++j) { ti += alpha[Row + j] * deltaP[j]; }
		quad += deltaP[i] * ti;     // δPᵀ α δP
		lin += deltaP[i] * beta[i]; // δPᵀ β
	}
	return 2.0 * lin - quad; // diff_pred
}
} // namespace palm_tracer
