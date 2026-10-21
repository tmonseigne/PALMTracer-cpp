#include "GaussianFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <limits>
#include <numbers>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "GaussianModel.hpp"
#include "signal_processing/Matrix.hpp"

namespace palm_tracer::fit {
namespace col = localization_columns;

namespace {
/// <summary>Prépare une fois les constantes du modèle courant ou candidat pour tous les pixels de la ROI.</summary>
GaussianModel PrepareModel(const fit_mode mode, const std::vector<double>& p, const double sigma, const double theta)
{
	switch (mode) {
		case fit_mode::gaussian_2d_fixed: return {p[0], p[1], sigma, false};
		case fit_mode::gaussian_2d_sigma: return {p[0], p[1], p[4], true};
		case fit_mode::gaussian_3d_fixed: return {p[0], p[1], p[4], p[5], theta, false};
		case fit_mode::gaussian_3d_theta: return {p[0], p[1], p[4], p[5], p[6], true};
		case fit_mode::none:
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: throw std::invalid_argument("GaussianFit : mode non gaussien.");
	}
	return {0, 0, 1, false}; // Pour éviter un warning de compilation, jamais atteint.
}
} // namespace

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
void GaussianFit::Fit(const localization::SegmentationInfo& infos, std::array<double, col::Count>& point) const
{
	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const auto Ux = static_cast<size_t>(infos.x_), Uy = static_cast<size_t>(infos.y_);
	// On est sur les bords, on ne fait pas de fitting, on remplit avec des valeurs par défaut (offset, intensity, sigmaX/Y, theta = -1)
	if (roiHalf_ > Ux || Ux >= dim_.x - roiHalf_ || roiHalf_ > Uy || Uy >= dim_.y - roiHalf_) {
		fill(infos, point, infos.x_, infos.y_, infos.z_);
		return;
	}

	// Valeurs par défaut
	const double Z = infos.z_;
	double x       = infos.x_, y     = infos.y_,
	       sigmaX  = sigmaX_, sigmaY = sigmaY_, theta = theta_;

	// Création des variables
	double lambda = LAMBDA_0, minI = 0, maxI = 0;
	// Un cache par thread couvre aussi les appels directs et les objets successifs de configurations différentes.
	static thread_local LM_Workspace ws;
	ws.Resize(n_, nParams_, nParams2_);

	// Initialisation des données (Récupération des min/max, de la ROI et initialisation des paramètres du modèle et de ses bornes)
	preCompute(Ux, Uy, ws.grid, minI, maxI);
	initP(ws.p, x, y, Z, minI, maxI);

	// Première itération de Levenberg-Marquardt et initialisation de l'erreur de référence
	double err = lm(ws);
	if (!std::isfinite(err) || !std::ranges::all_of(ws.beta, [](const double value) { return std::isfinite(value); })) {
		fill(infos, point, infos.x_, infos.y_, infos.z_);
		return;
	}

	// Trouver la valeur absolue maximale
	double m = std::ranges::max(std::views::transform(ws.beta, [](const double b) { return std::abs(b); }));

	for (size_t iteration = ((m < EPSILON_1) ? N_ITER : 0); iteration < N_ITER; ++iteration) { // Main Loop
		// Calcul de l'inverse de la matrice de courbure α
		ws.alphaR = ws.alpha; // Copie de α
		for (size_t i = 0; i < nParams_; ++i) {
			const size_t Diagonal = (i * nParams_) + i;
			// Une colonne de Jacobien nulle (angle à l'isotropie, par exemple) impose un pas nul.
			// Ce paramètre est réévalué à chaque itération et redevient libre dès qu'il est identifiable.
			ws.alphaR[Diagonal] = ws.alpha[Diagonal] == 0.0 ? 1.0 : (1.0 + lambda) * ws.alpha[Diagonal];
		}
		try { matrix::InverseGaussJordan(ws.alphaR, ws.alphaRi, nParams_); }
		catch (const std::domain_error&) {
			// Une régularisation plus forte peut rétablir un système numériquement inversible.
			if (lambda < LAMBDA_MAX) {
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX);
				continue;
			}
			fill(infos, point, infos.x_, infos.y_, infos.z_);
			return;
		}

		// Calcul du candidat (p + δp)
		std::ranges::fill(ws.deltaP, 0.0);
		for (size_t i = 0; i < nParams_; ++i) {
			for (size_t j = 0; j < nParams_; ++j) { ws.deltaP[i] += ws.alphaRi[(i * nParams_) + j] * ws.beta[j]; } // δP = αRi × β
			if (!std::isfinite(ws.deltaP[i])) {
				fill(infos, point, infos.x_, infos.y_, infos.z_);
				return;
			}
			ws.pTry[i] = std::clamp(ws.p[1][i] + ws.deltaP[i], ws.p[0][i], ws.p[2][i]); // pTry = clamp(p + δP)
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
			if (m < EPSILON_2 && iteration >= 2) { break; } // On est trop faible est déjà 2 itérations, on arrête
			continue;                                       // On repart au début de la boucle
		}

		// m = δPᵀ × (λ × δP + β)
		for (size_t i = 0; i < nParams_; ++i) { ws.lBetaDeltaP[i] = (lambda * ws.deltaP[i]) + ws.beta[i]; }
		m = std::inner_product(ws.deltaP.begin(), ws.deltaP.end(), ws.lBetaDeltaP.begin(), 0.0);

		// Calcul de l'erreur pour l'essai en cours
		const double Err_Try = getError(ws);
		const double Rho     = (err - Err_Try) / (2.0 * m);

		if (Rho > EPSILON_4) {                                  // On accepte : on remplace l'état courant par _try
			ws.p[1] = ws.pTry;                                  // Mise à jour des paramètres du modèle
			err     = lm(ws);                                   // Mise à jour de la ROI du modèle, la matrice alpha, le vecteur beta et de l'erreur.
			lambda  = std::max(lambda / LAMBDA_DN, LAMBDA_MIN); // Diminution de λ.
		}
		else {                                                 // On refuse : paramètres, coût et Jacobien courants restent cohérents.
			if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max pour λ, on arrête
			lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de λ.
		}

		// Critères d'arrêts
		m = 0.0;
		for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::abs(ws.deltaP[i]) / std::max(std::abs(ws.p[1][i]), 1e-12)); }
		if (iteration >= 2) {                                         // Si on est déjà à 2 itérations
			if (m < EPSILON_2) { break; }                             // Si le pas est trop faible
			if (err / static_cast<double>(n_) < EPSILON_3) { break; } // Si l'erreur est trop faible
			if (std::ranges::max(std::views::transform(ws.beta, [](const double b) { return std::abs(b); })) < EPSILON_1) { break; }
		}
	}

	// Mise à jour des différents éléments
	const double Offset = ws.p[1][0], Intensity = ws.p[1][1];
	x                   = ws.p[1][2];
	y                   = ws.p[1][3];


	const double RSS = computeR2(ws); // Calcul du R² sur XY
	if (RSS < EPSILON) {
		// Une corrélation au carré trop faible invalide le fit ; les paramètres de sortie gardent les valeurs d'échec.
		fill(infos, point, x, y, infos.z_);
		return;
	}

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
	fill(infos, point, x, y, Z, Offset, Intensity, sigmaX, sigmaY, theta, RSS);
}

//==================== Levenberg-Marquardt ====================
//-------------------------------------------------------------------------------------------------
double GaussianFit::lm(LM_Workspace& ws) const
{
	std::ranges::fill(ws.alpha, 0.0);
	std::ranges::fill(ws.beta, 0.0);
	double err          = 0.0;
	const auto Model    = PrepareModel(mode_, ws.p[1], sigmaX_, theta_);
	const bool Circular = mode_ == fit_mode::gaussian_3d_theta && std::abs(ws.p[1][4] - ws.p[1][5]) <= std::numeric_limits<double>::epsilon();

	for (size_t i = 0; i < n_; ++i) {
		const auto X = static_cast<double>(ws.grid[i].x), Y = static_cast<double>(ws.grid[i].y);
		ws.hat[i]    = Model.Evaluate(X - ws.p[1][2], Y - ws.p[1][3], ws.grad); // PSF théorique au pixel (x, y) selon les paramètres p
		ws.delta[i]  = ws.grid[i].v - ws.hat[i];                                // Résidu (∂) = mesure - prédiction
		ws.deltaW[i] = ws.weight[i] * ws.delta[i];                              // Résidu pondéré (W × ∂)
		err += ws.delta[i] * ws.deltaW[i];                                      // mise à jour du Chi² (somme des carrés pondérée W × ∂ × ∂)
		// L'angle reste indéterminé à l'isotropie, y compris en présence d'un résidu d'arrondi.
		if (Circular) { ws.grad[6] = 0.0; }

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
	return (2.0 * lin) - quad; // 2 × δPᵀ × β - δPᵀ × α × δP
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::getError(const LM_Workspace& ws) const
{
	double res       = 0.0;
	const auto Model = PrepareModel(mode_, ws.pTry, sigmaX_, theta_);
	for (size_t i = 0; i < n_; ++i) {
		const double Dx = ws.grid[i].v - Model.Evaluate(static_cast<double>(ws.grid[i].x) - ws.pTry[2], static_cast<double>(ws.grid[i].y) - ws.pTry[3]);
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
} // namespace palm_tracer::fit
