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
	//--------------------------------------------------------------------------------------------------------------------------------------------
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

	[[maybe_unused]] auto printGrad = [&](const std::vector<double>& g) {
		std::cout << "\t grad [";
		for (size_t i = 0; i < nParams_; ++i) { std::cout << g[i] << " "; }
		std::cout << "]\n";
	};
	//--------------------------------------------------------------------------------------------------------------------------------------------

	// Valeurs de sortie par défaut
	double offset = 0, intensity = 0, rss      = -1,
	       x      = infos.x_, y  = infos.y_, z = initZ_;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const auto Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux < dim_.x - roiHalf_ && roiHalf_ <= Uy && Uy < dim_.y - roiHalf_) {
		// Création des variables
		double lambda = LAMBDA_0, minI = 0, maxI = 0;
		LM_Workspace ws;
		ws.Resize(n_, nParams_, nParams2_);

		// Initialisation des données (Récupération des min/max, de la ROI et initialisation du P et de ses bornes)
		preCompute(Ux, Uy, ws.grid, minI, maxI);
		matrix::CenterOfMass(ws.grid, roi_, x, y); // Calcul du centre de masse de la région pour x et y ATTENTION !! x,y ∈ [0,roi]²
		initP(ws.p, x, y, z, minI, maxI);          // Initialisation du/des P
		// Amplitude (et offset) de départ (et max) cohérente avec le patch courant masi chelou cette etape, je ne visualise pas trop l'interet et si c'est négatif ça casse tout.
		//for (size_t i = 1; i <= 2; ++i) { for (size_t j = 0; j <= 1; ++j) { p[i][j] *= coeffRatio_; } }
		ws.p[2][4] = static_cast<double>(splineDim_.z);                                                        // Z ∈ [0, Zmax] en index de plan (flottant)
		for (size_t i = 0; i < nParams_; ++i) { ws.p[1][i] = std::clamp(ws.p[1][i], ws.p[0][i], ws.p[2][i]); } // Clamp en cas de problème.

		// Première itération de Levenberg-Marquardt
		double err     = lm(ws);
		double err_old = err;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(ws.beta, [](const double b) { return std::fabs(b); }));

		std::cout << std::setprecision(6) << std::fixed;
		std::cout << "\t [min .. val .. max]\n";
		constexpr std::array PName{"B", "I", "X", "Y", "Z"};
		for (size_t i = 0; i < nParams_; ++i) { std::cout << "\t" << PName[i] << "= [" << ws.p[0][i] << " .. " << ws.p[1][i] << " .. " << ws.p[2][i] << "]\n"; }
		printAlpha(ws.alpha);
		printBeta(ws.beta);
		printGrad(ws.grad);
		//--------------------------------------------------------------------------------------------------------------------------------------------

		for (size_t iteration = ((m < EPSILON_1) ? N_ITER : 0); iteration < N_ITER; ++iteration) { // Main Loop
			//--------------------------------------------------------------------------------------------------------------------------------------------
			std::cout << "\t Iteration " << iteration << "\t lambda = " << lambda << "\n";
			//--------------------------------------------------------------------------------------------------------------------------------------------

			// Calcul de l'inverse de la matrice de courbure α
			ws.alphaR = ws.alpha;                                                                                         // Copie de α
			for (size_t i = 0; i < nParams_; ++i) { ws.alphaR[i * nParams_ + i] += lambda * ws.alpha[i * nParams_ + i]; } // αR = α + λ diag(α)
			matrix::InverseGaussJordan(ws.alphaR, ws.alphaRi, nParams_);                                                  // Inversion (Gauss-Jordan)

			// Calcul du candidat (p + δp)
			std::ranges::fill(ws.deltaP, 0.0);
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { ws.deltaP[i] += ws.alphaRi[i * nParams_ + j] * ws.beta[j]; } // δP = αRi * β
				ws.p_try[i] = std::clamp(ws.p[1][i] + ws.deltaP[i], ws.p[0][i], ws.p[2][i]);                         // p_try = clamp(p + δP)
			}

			// Gain prédit (pas d’éval d’image ici) m = δPᵀ * (λ δP + β)
			const double Gain = gain(ws);
			if (Gain <= 0.0) {                                     // clairement pas bon → on n’évalue pas, on augmente λ et on retente
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de lambda
				// critères d’arrêt “petit pas”
				m = 0.0;
				for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::fabs(ws.deltaP[i]) / std::max(std::fabs(ws.p[1][i]), 1e-12)); }
				if (m < EPSILON_2 && iteration >= 2) { break; } // On est trop faible est déjà 2 iterations, on arrête
				continue;                                       // On repart au début de la boucle
			}

			// m = δPᵀ * (λ δP + β)
			for (size_t i = 0; i < nParams_; ++i) { ws.lBetaDeltaP[i] = lambda * ws.deltaP[i] + ws.beta[i]; }
			m = std::inner_product(ws.deltaP.begin(), ws.deltaP.end(), ws.lBetaDeltaP.begin(), 0.0);

			// On calcule les nouvelles valeurs : évalue image + lm pour p_try
			const double Err_Try = getErr(ws);
			const double Rho     = (err - Err_Try) / (2.0 * m);

			//--------------------------------------------------------------------------------------------------------------------------------------------
			std::cout << "\tPoint (" << ws.p[1][2] << ", " << ws.p[1][3] << ", " << ws.p[1][4]
				<< ")\n\t\t Iteration " << iteration << "\t err = " << err << "\t err_try = " << Err_Try << "\t m = " << m
				<< "\t Rho = " << Rho << "\t lambda = " << lambda << "\n";
			//--------------------------------------------------------------------------------------------------------------------------------------------

			if (Rho > EPSILON_4) {  // On accepte : on remplace l’état courant par _try
				err_old = err;      // Mise à jour de l'erreur de reference
				ws.p[1] = ws.p_try; // Mise à jour du modèle
				err     = lm(ws);   // Mise à jour de la ROI du modèle, la matrice alpha, le vecteur beta et de l'erreur.
				//--------------------------------------------------------------------------------------------------------------------------------------------
				std::cout << "\t\t Rho > EPSILON_4 Update Alpha and beta " << iteration << "\t lambda = " << lambda << "\n";
				printAlpha(ws.alpha);
				printBeta(ws.beta);
				printGrad(ws.grad);
				//--------------------------------------------------------------------------------------------------------------------------------------------
				lambda = std::max(lambda / LAMBDA_DN, LAMBDA_MIN); // Diminution du lambda
			}
			else {                                                 // On refuse : on garde l’état courant
				err = err_old;                                     // On reprend l'errur de reference
				if (lambda >= LAMBDA_MAX) { break; }               // On est déjà au max pour lambda, on arrête
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_MAX); // Augmentation de lambda
			}

			// Critères d'arrêts
			m = 0.0;
			for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::abs(ws.deltaP[i]) / std::max(std::abs(ws.p[1][i]), 1e-12)); }
			if (iteration >= 2) {                                         // Si on est déjà à 2 iterations
				if (m < EPSILON_2) { break; }                             // Si le pas est trop faible
				if (err / static_cast<double>(n_) < EPSILON_3) { break; } // Si l'erreur est trop faible
				if (std::ranges::max(std::views::transform(ws.beta, [](const double b) { return std::fabs(b); })) < EPSILON_1) { break; }
			}
		}

		// Assignation
		//offset = p[1][0] / coeffRatio_; // Retour à une intensité "normale"
		//intensity = p[1][1] / coeffRatio_; // Retour à une intensité "normale"
		offset    = ws.p[1][0];
		intensity = ws.p[1][1];
		x         = ws.p[1][2] - static_cast<double>(roiHalf_) + infos.x_;
		y         = ws.p[1][3] - static_cast<double>(roiHalf_) + infos.y_;
		z         = (ws.p[1][4] - initZ_); // INDEX DE PLAN (flottant), pas de rescaling (deltaZ) mais recentrage (p[1][4] - initZ_) donc Z ∈ [-initZ_ ; initZ_]
		rss       = (ws.beta[2] + ws.beta[3]) / 2;
	}
	else { z = 0; }

	std::cout << "Final Point : " << offset << ", " << intensity << ", " << x << ", " << y << ", " << z << "\n";
	fill(infos, point, offset, intensity, x, y, z, 0, 0, 0, rss, rss);
}

//-------------------------------------------------------------------------------------------------
void SplineFit::delta3D(LM_Workspace& ws) const
{
	// Coordonnées fractionnaires (ROI locale)
	const double Xr = -1.0 * (ws.p[1][2] - static_cast<double>(roiHalf_) + 0.5), Yr = -1.0 * (ws.p[1][3] - static_cast<double>(roiHalf_) + 0.5);

	// Parties entières pour l’ancrage du patch
	const int Xi = static_cast<int>(std::floor(Xr)), Yi = static_cast<int>(std::floor(Yr)), Zi = static_cast<int>(std::floor(ws.p[1][4]));

	// Parties fractionnaires à l’intérieur de la cellule
	const double Xc = Xr - Xi, Yc = Yr - Yi, Zc = ws.p[1][4] - Zi;

	std::ranges::fill(ws.f, 0.0);
	for (auto& df : ws.df) { std::ranges::fill(df, 0.0); }

	double cz = 1.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIDE; ++i) {
		double cy = 1.0;
		for (size_t j = 0; j < SPLINE_PATCH_SIDE; ++j) {
			double cx = 1.0;
			for (size_t k = 0; k < SPLINE_PATCH_SIDE; ++k) {
				const size_t Id   = i * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k;
				const double Coef = cz * cy * cx;
				ws.f[Id]          = Coef;
				if (k > 0) { ws.df[0][Id - 1] += static_cast<double>(k) * Coef; }                  // x : −1
				if (j > 0) { ws.df[1][Id - SPLINE_PATCH_SIDE] += static_cast<double>(j) * Coef; }  // y : −4
				if (i > 0) { ws.df[2][Id - SPLINE_PATCH_PLANE] += static_cast<double>(i) * Coef; } // z : −16
				cx *= Xc;
			}
			cy *= Yc;
		}
		cz *= Zc;
	}
}

//-------------------------------------------------------------------------------------------------
double SplineFit::lm(LM_Workspace& ws) const
{
	std::ranges::fill(ws.alpha, 0.0);
	std::ranges::fill(ws.beta, 0.0);
	double err = 0.0;

	delta3D(ws); // Mise à jour des contributions d'une base de fonctions splines 3D (valeurs et dérivées).
	// On se base sur les décimales donc f et df mis à jour ici ne changent pas pour toute la ROI
	// (ça à l'air d'etre une sorte de "poids" pour chacun des 64 coeffs)

	Vec3I pos{.x = 0, .y = 0, .z = static_cast<size_t>(ws.p[1][4])};
	for (size_t y = 0, idx = 0; y < roi_; ++y) {
		pos.y = static_cast<int>(y) + static_cast<int>(splineOffset_);
		for (size_t x = 0; x < roi_; ++x, ++idx) {
			pos.x          = static_cast<int>(x) + static_cast<int>(splineOffset_);
			ws.hat[idx]    = psf(ws, ws.p[1], pos);          // PSF théorique au pixel (x, y) selon les paramètres p
			ws.delta[idx]  = ws.grid[idx].v - ws.hat[idx];   // Résidu (∂) = mesure - prédiction
			ws.deltaW[idx] = ws.weight[idx] * ws.delta[idx]; // Résidu pondéré (W * ∂)
			err += ws.delta[idx] * ws.deltaW[idx];           // mise à jour du Chi² (somme des carrés pondérée W * ∂ * ∂)
			deriv(ws, pos);                                  // dérivées au pixel

			// Accumulations alpha += w * grad * gradᵀ  et  beta += w * grad * r
			for (size_t a = 0; a < nParams_; ++a) {
				ws.beta[a] += ws.grad[a] * ws.deltaW[idx];
				// rang-1 symétrique
				const size_t Row = a * nParams_;
				for (size_t b = 0; b < nParams_; ++b) { ws.alpha[Row + b] += ws.weight[idx] * ws.grad[a] * ws.grad[b]; }
			}
		}
	}
	return err;
}

//-------------------------------------------------------------------------------------------------
double SplineFit::gain(const LM_Workspace& ws) const
{
	// t = alpha * deltaP
	double quad = 0.0, lin = 0.0;
	for (size_t i = 0; i < nParams_; ++i) {
		double ti        = 0.0;
		const size_t Row = i * nParams_;
		for (size_t j = 0; j < nParams_; ++j) { ti += ws.alpha[Row + j] * ws.deltaP[j]; }
		quad += ws.deltaP[i] * ti;        // δPᵀ α δP
		lin += ws.deltaP[i] * ws.beta[i]; // δPᵀ β
	}
	return 2.0 * lin - quad; // diff_pred
}

//-------------------------------------------------------------------------------------------------
void SplineFit::deriv(LM_Workspace& ws, Vec3I& pos) const
{
	pos = pos.Clamp(splineDim_ - 1);

	const size_t Pos = pos.GetIndex(splineDim_);
	std::ranges::fill(ws.grad, 0.0); // Initialise toutes les dérivées à zéro

	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) {
		const size_t Id = i * splineSize_ + Pos;
		ws.grad[1] += ws.f[i] * coeffs_[Id];     // Contribution à l’intensité.
		ws.grad[2] += ws.df[0][i] * coeffs_[Id]; // Dérivée en x.
		ws.grad[3] += ws.df[1][i] * coeffs_[Id]; // Dérivée en y.
		ws.grad[4] += ws.df[2][i] * coeffs_[Id]; // Dérivée en z.
	}
	ws.grad[0] = 1;                // Dérivée par rapport à l’offset
	ws.grad[2] *= -1 * ws.p[1][1]; // Ajustement en x.
	ws.grad[3] *= -1 * ws.p[1][1]; // Ajustement en y.
	ws.grad[4] *= ws.p[1][1];      // Ajustement en z.
}

//-------------------------------------------------------------------------------------------------
double SplineFit::psf(const LM_Workspace& ws, const std::vector<double>& p, const Vec3I& pos) const
{
	const size_t Pos = pos.GetIndex(splineDim_);
	double s         = 0.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) { s += ws.f[i] * coeffs_[i * splineSize_ + Pos]; }
	return p[0] + p[1] * s;
}

//-------------------------------------------------------------------------------------------------
double SplineFit::getErr(const LM_Workspace& ws) const
{
	double res = 0.0;
	Vec3I pos{.x = 0, .y = 0, .z = static_cast<size_t>(ws.p_try[4])};
	for (size_t y = 0, idx = 0; y < roi_; ++y) {
		pos.y = static_cast<int>(y) + static_cast<int>(splineOffset_);
		for (size_t x = 0; x < roi_; ++x, ++idx) {
			pos.x           = static_cast<int>(x) + static_cast<int>(splineOffset_);
			const double Dx = ws.grid[idx].v - psf(ws, ws.p_try, pos); // PSF théorique au pixel (x, y) selon les paramètres p
			res += ws.weight[idx] * Dx * Dx;
		}
	}
	return res;
}
} // namespace palm_tracer
