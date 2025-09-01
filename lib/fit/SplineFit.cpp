#include "SplineFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <execution>
#include <iomanip>
#include <limits>
#include <ranges>
#include <vector>

#include "signal_processing/Matrix.hpp"

namespace palm_tracer {
/// <summary>Calcule la racine carrée réelle d'un nombre, même s'il est négatif.
///
/// Cette fonction retourne la partie réelle de la racine carrée complexe d'un nombre réel.
/// Contrairement à la fonction standard (std::sqrt(x)), qui retourne NaN pour une entrée négative.
/// Cette version utilise les nombres complexes pour permettre un traitement sûr.
///
/// Cela garantit :
/// - SqrtReal(x) = sqrt(x) si x ≥ 0
/// - SqrtReal(x) = 0 si x < 0 (car partie réelle de sqrt(x + 0i) est nulle pour x < 0)
/// </summary>
/// <param name="value">Nombre réel dont on veut extraire la racine carrée réelle.</param>
/// <returns>La partie réelle de la racine carrée complexe de `value`.</returns>
static double SqrtReal(const double value)
{
	const std::complex Z(value, 0.0);
	return std::real(std::sqrt(Z));
}

//-------------------------------------------------------------------------------------------------
void SplineFit::parseParams()
{
	size_t idx = 0;
	roi_       = static_cast<size_t>(fitParams_[idx++]);
	for (size_t i = 0; i < 3; ++i) { splineDim_[i] = static_cast<size_t>(fitParams_[idx++]); }
	splineSizeXY_ = splineDim_.GetPlaneSize();               // Taille d'un plan de coefficient : [Rx x Ry x Rz]
	splineSize_   = splineDim_.GetSize();                    // Taille du volume de coefficient : [Rx x Ry x Rz]
	initZ_        = static_cast<double>(splineDim_.z) / 2.0; // Centre de pile (utile pour l'initialisation)
	deltaZ_       = fitParams_[idx++];                       // Conservé pour compatibilité, pas utilisé dans la sortie
	coeffs_.resize(splineSize_ * SPLINE_PATCH_SIZE);         // Coefficients: [Rx x Ry x Rz] x 64 (4x4x4)
	for (double& coeff : coeffs_) { coeff = fitParams_[idx++]; }

	// Ratio d'intensité pour initialisation (normalisé sur le coeff central de la spline)
	const double CenterCoeff = coeffs_[32 * splineSize_ + (splineDim_.z / 2) * splineSizeXY_ + (splineDim_.y / 2) * splineDim_.x + (splineDim_.x / 2)];
	std::cout << "Center Coeff = " << CenterCoeff << "\n";
	coeffRatio_ = std::fabs(CenterCoeff) < std::numeric_limits<double>::epsilon() ? 1.0 : 4.0 / CenterCoeff;
	// Décalage pour recaler la ROI sur la grille de spline (en cas de tailles différentes)
	splineOffset_ = static_cast<size_t>(std::floor((static_cast<double>(splineDim_.x) + 1.0 - static_cast<double>(roi_)) / 2));
}

//-------------------------------------------------------------------------------------------------
void SplineFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const
{
	[[maybe_unused]] auto printBeta = [&](const std::vector<double>& b) {
		std::cout << "\t beta [";
		for (size_t i = 0; i < nParams_; ++i) { std::cout << b[i] << " "; }
		std::cout << "]\n";
	};

	[[maybe_unused]] auto printAlpha = [&](const std::vector<double>& a) {
		std::cout << "\t Alpha : \n";
		for (size_t i = 0, id = 0; i < nParams_; ++i) {
			std::cout << "\t| ";
			for (size_t j = 0; j < nParams_; ++j, ++id) { std::cout << a[id] << " | "; }
			std::cout << "\n";
		}
		std::cout << "\n";
	};


	// Valeurs de sortie par défaut
	double offset = 0, intensity = 0,
	       x      = infos.x_, y  = infos.y_, z = initZ_,
	       rssXY  = -1, rssZ     = -1;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const auto Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux < dim_.x - roiHalf_ && roiHalf_ <= Uy && Uy < dim_.y - roiHalf_) {
		// Création des variables
		double lambda = LAMBDA_0, minI = 0, maxI = 0;
		std::vector<double> alpha(nParams2_), alphaR(nParams2_), alphaRi(nParams2_),
		                    beta(nParams_), deltaP(nParams_), lBetaDeltaP(nParams_), p_try(nParams_), grad(nParams_);
		std::vector<Pixel> grid(n_);
		std::array<std::vector<double>, 3> p;
		std::vector<double> f(SPLINE_PATCH_SIZE);
		std::array df{std::vector<double>(SPLINE_PATCH_SIZE), std::vector<double>(SPLINE_PATCH_SIZE), std::vector<double>(SPLINE_PATCH_SIZE)};

		// Initialisation des données (Récupération des min/max, de la ROI et initialisation du P et de ses bornes)
		preCompute(Ux, Uy, grid, minI, maxI);
		matrix::CenterOfMass(grid, roi_, x, y); // Calcul du centre de masse de la région pour x et y
		initP(p, x, y, z, minI, maxI);          // Initialisation du/des P
		// Amplitude (et offset) de départ (et max) cohérente avec le patch courant masi chelou cette etape, je ne visualise pas trop l'interet et si c'est négatif ça casse tout.
		//for (size_t i = 1; i <= 2; ++i) { for (size_t j = 0; j <= 1; ++j) { p[i][j] *= coeffRatio_; } }
		p[2][4] = static_cast<double>(splineDim_.z);                                               // Z ∈ [0, Zmax] en index de plan (flottant)
		for (size_t i = 0; i < nParams_; ++i) { p[1][i] = std::clamp(p[1][i], p[0][i], p[2][i]); } // Clamp en cas de problème.

		// Première itération de Levenberg-Marquardt
		double err     = lm(grid, p[1], alpha, beta, f, df, grad, 0);
		double err_old = err;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

		std::cout << std::setprecision(6) << std::fixed;
		std::cout << "\t [min .. val .. max]\n";
		constexpr std::array PName{"B", "I", "X", "Y", "Z"};
		for (size_t i = 0; i < nParams_; ++i) { std::cout << "\t" << PName[i] << "= [" << p[0][i] << " .. " << p[1][i] << " .. " << p[2][i] << "]\n"; }
		printAlpha(alpha);
		printBeta(beta);
		//--------------------------------------------------------------------------------------------------------------------------------------------

		const bool Stop = (m < EPSILON_1);

		for (size_t iteration = (Stop ? N_ITER : 0); iteration < N_ITER; ++iteration) { // Main Loop
			alphaR = alpha;
			//for (size_t i = 0; i < nParams_; ++i) { alphaR[i * nParams_ + i] += lambda; } // alphaR = alpha + λ (amortissement Marquardt, robuste si diag ~ 0)
			for (size_t i = 0; i < nParams_; ++i) { alphaR[i * nParams_ + i] += lambda * alpha[i * nParams_ + i]; } // alphaR = alpha + λ diag(alpha)
			matrix::InverseGaussJordan(alphaR, alphaRi, nParams_);                                                  // Inversion (Gauss-Jordan) et pas LM

			std::ranges::fill(deltaP, 0.0);
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { deltaP[i] += alphaRi[i * nParams_ + j] * beta[j]; } // δP
				p_try[i] = std::clamp(p[1][i] + deltaP[i], p[0][i], p[2][i]);                               // p_try = clamp(p + δP)
			}

			// Gain prédit (deltaPᵀ * (λ deltaP + beta)) et diff prédite (2 lin − quad)
			for (size_t i = 0; i < nParams_; ++i) { lBetaDeltaP[i] = lambda * deltaP[i] + beta[i]; }
			m = std::inner_product(deltaP.begin(), deltaP.end(), lBetaDeltaP.begin(), 0.0);
			if (!std::isfinite(m) || m <= 0.0) {
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX);
				continue; // on réessaie avec plus de damping
			}


			const double Gain = gain(deltaP, alpha, beta);
			if (Gain <= 0.0) {                                     // clairement pas bon → on n’évalue pas, on augmente λ et on retente
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de lambda

				// critères d’arrêt “petit pas”
				m = 0.0;
				for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::fabs(deltaP[i]) / std::max(std::fabs(p[1][i]), 1e-12)); }
				if (m < EPSILON_2 && iteration >= 2) { break; } // On est trop faible est déjà 2 iterations, on arrête
				continue;                                       // On repart au début de la boucle
			}

			// Évalue l'erreur au point d'essai p_try (sans recalculer alpha/beta)
			const double Err_Try = lm(grid, p_try, alpha, beta, f, df, grad, 1);
			if (!std::isfinite(Err_Try)) {
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX);
				continue;
			}

			const double Rho = (err - Err_Try) / (2.0 * m);

			std::cout << "\tPoint (" << p[1][2] << ", " << p[1][3] << ", " << p[1][4]
				<< ")\n\t\t Iteration " << iteration << "\t err = " << err << "\t err_try = " << Err_Try << "\t m = " << m
				<< "\t Rho = " << Rho << "\t lambda = " << lambda << "\n";

			if (Rho > EPSILON_4) { // On accepte : on remplace l'état courant et on recalcule alpha/beta
				err_old = err;
				p[1]    = p_try;
				err     = lm(grid, p[1], alpha, beta, f, df, grad, 0);
				std::cout << "\t\tUpdate: err = " << err << "\n";
				printAlpha(alpha);
				printBeta(beta);
				lambda = std::max(lambda / LAMBDA_DN, LAMBDA_MIN);
			}
			else { // On refuse : on garde l'état courant et on augmente λ
				err = err_old;
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de lambda
			}

			// Critères d'arrêt après quelques itérations
			m = 0.0;
			for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::fabs(deltaP[i]) / std::max(std::fabs(p[1][i]), 1e-12)); }
			if (iteration >= 2) { // Si on est déjà à 2 iterations, on arrête
				if (m < EPSILON_2) { break; }
				if (err / static_cast<double>(n_) < EPSILON_3) { break; }
				if (std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); })) < EPSILON_1) { break; }
			}
		}

		lm(grid, p[1], alpha, beta, f, df, grad, 2);

		// Assignation
		//offset = p[1][0] / coeffRatio_; // Retour à une intensité "normale"
		//intensity = p[1][1] / coeffRatio_; // Retour à une intensité "normale"
		offset    = p[1][0];
		intensity = p[1][1];
		x         = p[1][2] - static_cast<double>(roiHalf_) + infos.x_;
		y         = p[1][3] - static_cast<double>(roiHalf_) + infos.y_;
		z         = (p[1][4] - initZ_); // INDEX DE PLAN (flottant), pas de rescaling (deltaZ) mais recentrage (p[1][4] - initZ_) donc Z ∈ [-initZ_ ; initZ_]
		rssXY     = (beta[2] + beta[3]) / 2;
		rssZ      = beta[4];
	}
	else { z = 0; }

	std::cout << "Final Point : " << offset << ", " << intensity << ", " << x << ", " << y << ", " << z << "\n";
	fill(infos, point, offset, intensity, x, y, z, 0, 0, 0, rssXY, rssZ);
}

//-------------------------------------------------------------------------------------------------
void SplineFit::derivativeSpline(Vec3I& pos, const std::vector<double>& f, const std::array<std::vector<double>, 3>& df,
                                 const std::vector<double>& p, std::vector<double>& dudt, double& model) const
{
	pos = pos.Clamp(splineDim_ - 1);

	const size_t Pos = pos.GetIndex(splineDim_);
	std::ranges::fill(dudt, 0.0); // Initialise toutes les dérivées à zéro


	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) {
		const size_t Id = i * splineSize_ + Pos;
		dudt[1] += f[i] * coeffs_[Id];     // Contribution à l’intensité.
		dudt[2] += df[0][i] * coeffs_[Id]; // Dérivée en x.
		dudt[3] += df[1][i] * coeffs_[Id]; // Dérivée en y.
		dudt[4] += df[2][i] * coeffs_[Id]; // Dérivée en z.
	}
	dudt[0] = 1;                   // Dérivée par rapport à l’offset
	dudt[2] *= -1 * p[1];          // Ajustement en x.
	dudt[3] *= -1 * p[1];          // Ajustement en y.
	dudt[4] *= p[1];               // Ajustement en z.
	model = p[0] + p[1] * dudt[1]; // Évaluation du modèle spline : offset + amplitude * contribution
}

//-------------------------------------------------------------------------------------------------
void SplineFit::delta3D(const std::vector<double>& p, std::vector<double>& f, std::array<std::vector<double>, 3>& df, int& xi, int& yi, int& zi) const
{
	const double Xr = -1.0 * (p[2] - static_cast<double>(roiHalf_) + 0.5), Yr = -1.0 * (p[3] - static_cast<double>(roiHalf_) + 0.5);

	xi = static_cast<int>(std::floor(Xr));
	yi = static_cast<int>(std::floor(Yr));
	zi = static_cast<int>(std::floor(p[4]));

	const double Xc = Xr - xi, Yc = Yr - yi, Zc = p[4] - zi;

	std::ranges::fill(f, 0.0);
	for (size_t i = 0; i < 3; ++i) { std::ranges::fill(df[i], 0.0); }
	double cz = 1.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIDE; ++i) {
		double cy = 1.0;
		for (size_t j = 0; j < SPLINE_PATCH_SIDE; ++j) {
			double cx = 1.0;
			for (size_t k = 0; k < SPLINE_PATCH_SIDE; ++k) {
				const size_t Id   = i * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k;
				const double Coef = cz * cy * cx;
				f[Id]             = Coef;
				if (k > 0) { df[0][Id - 1] += static_cast<double>(k) * Coef; }                  // x : −1
				if (j > 0) { df[1][Id - SPLINE_PATCH_SIDE] += static_cast<double>(j) * Coef; }  // y : −4
				if (i > 0) { df[2][Id - SPLINE_PATCH_PLANE] += static_cast<double>(i) * Coef; } // z : −16
				cx *= Xc;
			}
			cy *= Yc;
		}
		cz *= Zc;
	}
}

//-------------------------------------------------------------------------------------------------
double SplineFit::lm(const std::vector<Pixel>& grid, const std::vector<double>& p, std::vector<double>& alpha, std::vector<double>& beta,
                     std::vector<double>& f, std::array<std::vector<double>, 3>& df, std::vector<double>& grad, const int mode) const
{
	std::ranges::fill(alpha, 0.0);
	std::ranges::fill(beta, 0.0);
	double res = 0.0;

	// Calcul du modèle pour chaque pixel de la ROI selon les paramètres p
	// Résidu (∂) = mesure - prédiction
	// Résidu pondéré (W * ∂)
	// mise à jour du Chi² (somme des carrés pondérée W * ∂ * ∂)
	// dérivées au pixel
	// Accumulations alpha += w * grad * gradᵀ  et  beta += w * grad * r

	std::vector mInv(mode == 2 ? nParams2_ : 0, 0.0);

	int xi, yi, zi;
	delta3D(p, f, df, xi, yi, zi);

	Vec3I pos{.x = 0, .y = 0, .z = static_cast<size_t>(zi)};
	for (size_t i = 0, idx = 0; i < roi_; ++i) {
		pos.x = xi + static_cast<int>(i) + static_cast<int>(splineOffset_);
		for (size_t j = 0; j < roi_; ++j, ++idx) {
			pos.y = yi + static_cast<int>(j) + static_cast<int>(splineOffset_);

			double model = 0.0;
			std::ranges::fill(grad, 0.0);
			derivativeSpline(pos, f, df, p, grad, model); //Calcul les dérivées du modèle spline 3D en un point donné (x, y, z).
			const double Mu = std::max(model, EPSILON_2);

			if (mode == 2) { // fisher
				// Matrice d'information de Fisher (approx)
				for (size_t a = 0; a < nParams_; ++a) {
					for (size_t b = a; b < nParams_; ++b) {
						alpha[a * nParams_ + b] += grad[b] * grad[a] / Mu;
						alpha[b * nParams_ + a] = alpha[a * nParams_ + b];
					}
				}
			}
			else {
				// Calcul du terme d’erreur : -2 * log-vraisemblance
				const double Data = grid[idx].v;
				if (Data > 0) { res += 2 * (Mu - Data - Data * std::log(Mu / Data)); }
				else { res += 2 * Mu; }
				// Mise à jour erreur/gradient/hessienne (EMCCD)
				if (mode == 0) {
					const double T1 = 1 - Data / Mu;
					for (size_t k = 0; k < nParams_; ++k) { beta[k] += T1 * grad[k]; } // ∂L/∂θ_i

					const double T2 = Data / (Mu * Mu);
					for (size_t a = 0; a < nParams_; ++a) {
						for (size_t b = a; b < nParams_; ++b) {
							alpha[a * nParams_ + b] += T2 * grad[a] * grad[b];
							alpha[b * nParams_ + a] = alpha[a * nParams_ + b]; // Hessienne symétrique
						}
					}
				}
			}
		}
	}

	if (mode == 2) {
		for (size_t k = 0; k < nParams_; ++k) {
			if (std::fabs(alpha[k * nParams_ + k]) < std::numeric_limits<double>::epsilon()) { alpha[k * nParams_ + k] = EPSILON_2; }
		}

		// Matrice inverse (CRLB=F⁻¹)
		matrix::InverseLU(alpha, mInv, nParams_);

		for (size_t k = 0; k < nParams_; ++k) { beta[k] = SqrtReal(mInv[k * nParams_ + k]); }
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
double SplineFit::gain(const std::vector<double>& deltaP, const std::vector<double>& alpha, const std::vector<double>& beta) const
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
