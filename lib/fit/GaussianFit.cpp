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
//==================== Common Functions ====================
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
		LM_Workspace ws;
		ws.Resize(n_, nParams_, nParams2_);

		// Initialisation des données (Récupération des min/max, de la ROI et initialisation des paramètres du modèle et de ses bornes)
		preCompute(Ux, Uy, ws.grid, minI, maxI);
		initP(ws.p, x, y, z, minI, maxI);

		// Première itération de Levenberg-Marquardt et initialisation de l'erreur de référence
		double err     = lm(ws);
		double err_old = err;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(ws.beta, [](const double b) { return std::abs(b); }));

		for (size_t iteration = ((m < EPSILON_1) ? N_ITER : 0); iteration < N_ITER; ++iteration) { // Main Loop
			// Calcul de l'inverse de la matrice de courbure α
			ws.alphaR = ws.alpha;                                                                                         // Copie de α
			for (size_t i = 0; i < nParams_; ++i) { ws.alphaR[i * nParams_ + i] += lambda * ws.alpha[i * nParams_ + i]; } // αR = α + λ × diag(α)
			matrix::InverseGaussJordan(ws.alphaR, ws.alphaRi, nParams_);                                                  // Inversion (Gauss-Jordan)

			// Calcul du candidat (p + δp)
			std::ranges::fill(ws.deltaP, 0.0);
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { ws.deltaP[i] += ws.alphaRi[i * nParams_ + j] * ws.beta[j]; } // δP = αRi × β
				ws.pTry[i] = std::clamp(ws.p[1][i] + ws.deltaP[i], ws.p[0][i], ws.p[2][i]);                          // pTry = clamp(p + δP)
			}
			// En cas de clamp de pTry, deltaP doit être recalculé.
			for (size_t i = 0; i < nParams_; ++i) { ws.deltaP[i] = ws.pTry[i] - ws.p[1][i]; }

			// Critère d'arrêt “petit pas”, si le pas réel est ~0 pour tous les paramètres. Inutile d'évaluer, on augmente λ.
			bool tinyStep = true;
			for (const double V : ws.deltaP) { tinyStep &= (std::abs(V) < 1e-15); }
			if (tinyStep) {
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de λ.
				continue;                                          // On repart au début de la boucle
			}

			// Gain prédit (pas d'éval d'image ici)
			const double Gain = gain(ws);
			if (Gain <= 0.0) {                                     // clairement pas bon → on n'évalue pas, on augmente λ et on retente
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de λ.
				// Critère d'arrêt “petit pas”, m = max (δP / p) (variation relative maximum parmi tous les paramètres)
				m = 0.0;
				for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::abs(ws.deltaP[i]) / std::max(std::abs(ws.p[1][i]), 1e-12)); }
				if (m < EPSILON_2 && iteration >= 2) { break; } // On est trop faible est déjà 2 iterations, on arrête
				continue;                                       // On repart au début de la boucle
			}

			// m = δPᵀ × (λ × δP + β)
			for (size_t i = 0; i < nParams_; ++i) { ws.lBetaDeltaP[i] = lambda * ws.deltaP[i] + ws.beta[i]; }
			m = std::inner_product(ws.deltaP.begin(), ws.deltaP.end(), ws.lBetaDeltaP.begin(), 0.0);

			// Calcul de l'erreur pour l'essai en cours
			const double Err_Try = getError(ws);
			const double Rho     = (err - Err_Try) / (2.0 * m);

			if (Rho > EPSILON_4) {                                  // On accepte : on remplace l'état courant par _try
				err_old = err;                                      // Mise à jour de l'erreur de reference
				ws.p[1] = ws.pTry;                                  // Mise à jour des paramètres du modèle
				err     = lm(ws);                                   // Mise à jour de la ROI du modèle, la matrice alpha, le vecteur beta et de l'erreur.
				lambda  = std::max(lambda / LAMBDA_DN, LAMBDA_MIN); // Diminution de λ.
			}
			else {                                                 // On refuse : on garde l'état courant
				err = err_old;                                     // On reprend l'erreur de reference
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max pour λ, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de λ.
			}

			// Critères d'arrêts
			m = 0.0;
			for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::abs(ws.deltaP[i]) / std::max(std::abs(ws.p[1][i]), 1e-12)); }
			if (iteration >= 2) {                                         // Si on est déjà à 2 iterations
				if (m < EPSILON_2) { break; }                             // Si le pas est trop faible
				if (err / static_cast<double>(n_) < EPSILON_3) { break; } // Si l'erreur est trop faible
				if (std::ranges::max(std::views::transform(ws.beta, [](const double b) { return std::abs(b); })) < EPSILON_1) { break; }
			}
		}

		// Mise à jour des différents éléments
		offset    = ws.p[1][0];
		intensity = ws.p[1][1];
		x         = ws.p[1][2];
		y         = ws.p[1][3];
		switch (mode_) {
			case fit_mode::none:
			case fit_mode::spline_3d_emccd:
			case fit_mode::spline_3d_scmos:
			case fit_mode::gaussian_2d_fixed: break;
			case fit_mode::gaussian_2d_sigma: {
				sigmaX = ws.p[1][4];
				break;
			}
			case fit_mode::gaussian_3d_fixed: {
				sigmaX = ws.p[1][4] / std::numbers::sqrt2_v<double>;
				sigmaY = ws.p[1][5] / std::numbers::sqrt2_v<double>;
				break;
			}
			case fit_mode::gaussian_3d_theta: {
				sigmaX = ws.p[1][4] / std::numbers::sqrt2_v<double>;
				sigmaY = ws.p[1][5] / std::numbers::sqrt2_v<double>;
				theta  = ws.p[1][6];
				break;
			}
		}
		rss = computeR2(ws); // Calcul du R² sur XY
	}

	// Post-normalisation selon le mode.
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

//==================== Levenberg-Marquardt ====================
//-------------------------------------------------------------------------------------------------
double GaussianFit::lm(LM_Workspace& ws) const
{
	std::ranges::fill(ws.alpha, 0.0);
	std::ranges::fill(ws.beta, 0.0);
	double err = 0.0;

	for (size_t i = 0; i < n_; ++i) {
		const double X = static_cast<double>(ws.grid[i].x), Y = static_cast<double>(ws.grid[i].y);
		ws.hat[i]      = psf(ws.p[1], X, Y);         // PSF théorique au pixel (x, y) selon les paramètres p
		ws.delta[i]    = ws.grid[i].v - ws.hat[i];   // Résidu (∂) = mesure - prédiction
		ws.deltaW[i]   = ws.weight[i] * ws.delta[i]; // Résidu pondéré (W × ∂)
		err += ws.delta[i] * ws.deltaW[i];           // mise à jour du Chi² (somme des carrés pondérée W × ∂ × ∂)
		deriv(ws.grad, ws.p[1], X, Y);               // dérivées au pixel

		// Accumulations alpha += w × grad × gradᵀ  et  beta += w × grad × r
		for (size_t a = 0; a < nParams_; ++a) {
			ws.beta[a] += ws.grad[a] * ws.deltaW[i];
			// rang-1 symétrique
			const size_t Row = a * nParams_;
			for (size_t b = 0; b < nParams_; ++b) { ws.alpha[Row + b] += ws.weight[i] * ws.grad[a] * ws.grad[b]; }
		}
	}
	return err;
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::gain(const LM_Workspace& ws) const
{
	double quad = 0.0, lin = 0.0;
	for (size_t i = 0; i < nParams_; ++i) {
		double ti        = 0.0;
		const size_t Row = i * nParams_;
		for (size_t j = 0; j < nParams_; ++j) { ti += ws.alpha[Row + j] * ws.deltaP[j]; }
		quad += ws.deltaP[i] * ti;        // δPᵀ × α × δP
		lin += ws.deltaP[i] * ws.beta[i]; // δPᵀ × β
	}
	return 2.0 * lin - quad; // 2 × δPᵀ × β - δPᵀ × α × δP
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::getError(const LM_Workspace& ws) const
{
	double res = 0.0;
	for (size_t i = 0; i < n_; ++i) {
		const double Dx = ws.grid[i].v - psf(ws.pTry, static_cast<double>(ws.grid[i].x), static_cast<double>(ws.grid[i].y));
		res += ws.weight[i] * Dx * Dx;
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::computeR2(LM_Workspace& ws) const
{
	// Calcul de la moyenne des données observées et du modèle
	const double MuD = std::accumulate(ws.grid.begin(), ws.grid.end(), 0.0, [](double sum, const auto& p) { return sum + p.v; }) / static_cast<double>(n_),
	             MuM = std::accumulate(ws.hat.begin(), ws.hat.end(), 0.0) / static_cast<double>(n_);

	// Normalisation en place (les éléments ne seront plus utilisés après)
	for (size_t i = 0; i < n_; ++i) {
		ws.grid[i].v -= MuD;
		ws.hat[i] -= MuM;
	}

	// Calcul de R = (D·M) / (||D|| × ||M||)
	const double T1 = std::inner_product(ws.grid.begin(), ws.grid.end(), ws.hat.begin(), 0.0, std::plus(),
	                                     [](const auto& a, double b) { return a.v * b; }), // D × M
	             T2 = std::inner_product(ws.grid.begin(), ws.grid.end(), ws.grid.begin(), 0.0, std::plus(),
	                                     [](const auto& a, const auto& b) { return a.v * b.v; }), // D × D
	             T3 = std::inner_product(ws.hat.begin(), ws.hat.end(), ws.hat.begin(), 0.0),      // M × M
	             R  = T1 / std::sqrt(T2 * T3);                                                    // Calcul du R

	return std::isnan(R) ? 0 : R * R;
}

// ==================== Run Computes ====================
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

//-------------------------------------------------------------------------------------------------
void GaussianFit::deriv(std::vector<double>& grad, const std::vector<double>& p, const double x, const double y) const
{
	switch (mode_) {
		case fit_mode::none:
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos:
		case fit_mode::gaussian_2d_fixed: Deriv(grad, p[1], x - p[2], y - p[3], sigmaX_, false);
			break;
		case fit_mode::gaussian_2d_sigma: Deriv(grad, p[1], x - p[2], y - p[3], p[4], true);
			break;
		case fit_mode::gaussian_3d_fixed: Deriv(grad, p[1], x - p[2], y - p[3], p[4], p[5], theta_, false);
			break;
		case fit_mode::gaussian_3d_theta: Deriv(grad, p[1], x - p[2], y - p[3], p[4], p[5], p[6], true);
			break;
	}
}
} // namespace palm_tracer
