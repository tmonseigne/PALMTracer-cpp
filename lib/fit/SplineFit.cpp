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
	coeffRatio_              = std::fabs(CenterCoeff) < std::numeric_limits<double>::epsilon() ? 1.0 : 4.0 / CenterCoeff;
	// Décalage pour recaler la ROI sur la grille de spline (en cas de tailles différentes)
	splineOffset_ = static_cast<size_t>(std::floor((static_cast<double>(splineDim_.x) + 1.0 - static_cast<double>(roi_)) / 2));

	jumpXY_ = 0.25;
	jumpZ_  = 0.25; // std::max(static_cast<double>(splineDim_.z) / 10.0, 0.5);
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
		std::array<double, 5> col{};
		for (size_t j = 0; j < nParams_; ++j) { col[j] = std::sqrt(std::abs(a[j * nParams_ + j])); }
		std::cout << "\t ||Alpha|| ≈ [" << col[0] << " " << col[1] << " " << col[2] << " " << col[3] << " " << col[4] << "]\n";
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
		ws.Imin = std::max(1.0, 0.01 * (maxI - minI)); // plancher intensité
		matrix::CenterOfMass(ws.grid, roi_, x, y);     // Calcul du centre de masse de la région pour x et y ATTENTION !! x,y ∈ [0,roi]²
		initP(ws.p, x, y, z, minI, maxI);              // Initialisation du/des P
		// Amplitude (et offset) de départ (et max) cohérente avec le patch courant masi chelou cette etape, je ne visualise pas trop l'interet et si c'est négatif ça casse tout.
		//for (size_t i = 1; i <= 2; ++i) { for (size_t j = 0; j <= 1; ++j) { p[i][j] *= coeffRatio_; } }
		ws.p[0][1] = std::max(ws.p[0][1], ws.Imin);                                                            // borne basse d'I
		ws.p[1][1] = std::max(ws.p[1][1], ws.Imin);                                                            // valeur courante ≥ I_min
		ws.p[2][4] = static_cast<double>(splineDim_.z);                                                        // Z ∈ [0, Zmax] en index de plan (flottant)
		for (size_t i = 0; i < nParams_; ++i) { ws.p[1][i] = std::clamp(ws.p[1][i], ws.p[0][i], ws.p[2][i]); } // Clamp en cas de problème.

		computeCoeffRatio(ws, ws.p[1]);

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
			ws.alphaR = ws.alpha; // Copie de α
			// αR = α + λ diag(α)
			for (size_t i = 0; i < nParams_; ++i) {
				const size_t I = i * nParams_ + i;
				const double D = ws.alpha[I];
				ws.alphaR[I] += lambda * ((std::fabs(D) > 1e-12) ? D : 1.0); // évite "λ * 0"
			}

			// scaling colonne (Marquardt "scaled")
			for (size_t j = 0; j < nParams_; ++j) { ws.sCol[j] = 1.0 / std::max(std::sqrt(std::abs(ws.alphaR[j * nParams_ + j])), 1e-12); }

			//αS = S * αR * S
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { ws.alphaS[i * nParams_ + j] = ws.sCol[i] * ws.alphaR[i * nParams_ + j] * ws.sCol[j]; }
			}

			//βS = S * -β
			for (size_t j = 0; j < nParams_; ++j) { ws.betaS[j] = -ws.sCol[j] * ws.beta[j]; }
			matrix::InverseGaussJordan(ws.alphaS, ws.alphaSi, nParams_); // Inversion (Gauss-Jordan)

			// Calcul du candidat (p + δp)
			bool clamped = false;
			std::ranges::fill(ws.deltaP, 0.0);
			for (size_t i = 0; i < nParams_; ++i) {
				for (size_t j = 0; j < nParams_; ++j) { ws.deltaP[i] += ws.alphaSi[i * nParams_ + j] * ws.betaS[j]; } // δP = αSi * βS
				ws.deltaP[i] *= ws.sCol[i];
				const double DPI = ws.deltaP[i];
				if (i == 2 || i == 3) { ws.deltaP[i] = std::clamp(ws.deltaP[i], -jumpXY_, jumpXY_); }
				else if (i == 4) { ws.deltaP[i] = std::clamp(ws.deltaP[i], -jumpZ_, jumpZ_); }
				else if (i == 1) {
					const double Jump = 0.5 * std::max(ws.p[1][1], ws.Imin); // pas max = ±50% de I (ou du plancher)
					ws.deltaP[i]      = std::clamp(ws.deltaP[i], -Jump, Jump);
				}
				if (std::fabs(DPI - ws.deltaP[i]) > 0.0) { clamped = true; }
				ws.p_try[i] = std::clamp(ws.p[1][i] + ws.deltaP[i], ws.p[0][i], ws.p[2][i]); // p_try = clamp(p + δP)
			}

			// Les splines partent vite, donc on prend le pas réel final (après clamp du p_try).
			for (size_t i = 0; i < nParams_; ++i) { ws.deltaP[i] = ws.p_try[i] - ws.p[1][i]; }
			// Petit garde-fou : si le pas réel est ~0, inutile d’évaluer
			bool tiny = true;
			for (const double V : ws.deltaP) { tiny &= (std::abs(V) < 1e-15); }
			if (tiny) {
				if (lambda >= LAMBDA_MAX) { break; }                                                              // On est déjà au max, on arrête
				lambda = std::min(lambda * (clamped ? LAMBDA_UP * std::sqrt(LAMBDA_UP) : LAMBDA_UP), LAMBDA_MAX); // Augmentation de lambda
				continue;                                                                                         // On repart au début de la boucle
			}

			// m = δPᵀ * (λ diag(α) δP + β) (diag(α)^pour retourner dnas l'échelle courante)
			for (size_t i = 0; i < nParams_; ++i) { ws.lBetaDeltaP[i] = lambda * ws.alpha[i * nParams_ + i] * ws.deltaP[i] - ws.beta[i]; }
			m = std::inner_product(ws.deltaP.begin(), ws.deltaP.end(), ws.lBetaDeltaP.begin(), 0.0);

			// Calcul des nouvelles valeurs : évalue image + lm pour p_try
			const double Err_Try = getErr(ws);
			const double Rho     = (err - Err_Try) / (2.0 * m);

			//--------------------------------------------------------------------------------------------------------------------------------------------
			std::cout << "\tPoint (" << ws.p[1][2] << ", " << ws.p[1][3] << ", " << ws.p[1][4]
				<< ")\n\t\t Iteration " << iteration << "\t err = " << err << "\t err_try = " << Err_Try << "\t m = " << m
				<< "\t Rho = " << Rho << "\t lambda = " << lambda << "\t coeffRatio = " << ws.coeffRatio << "\n";
			//--------------------------------------------------------------------------------------------------------------------------------------------

			if (err < Err_Try) { // REJET inconditionnel si l’erreur augmente (ou stagne)
				if (lambda >= LAMBDA_MAX) { break; }
				lambda = std::min(lambda * (clamped ? LAMBDA_UP * std::sqrt(LAMBDA_UP) : LAMBDA_UP), LAMBDA_MAX);
				continue;
			}

			if (m > 0.0 && Rho > EPSILON_4) { // On accepte : on remplace l’état courant par _try
				err_old = err;                // Mise à jour de l'erreur de reference
				ws.p[1] = ws.p_try;           // Mise à jour du modèle
				err     = lm(ws);             // Mise à jour de la ROI du modèle, la matrice alpha, le vecteur beta et de l'erreur.
				//--------------------------------------------------------------------------------------------------------------------------------------------
				std::cout << "\t\t Rho > EPSILON_4 Update Alpha and beta " << iteration << "\t lambda = " << lambda << "\n";
				printAlpha(ws.alpha);
				printBeta(ws.beta);
				printGrad(ws.grad);
				//--------------------------------------------------------------------------------------------------------------------------------------------
				lambda = std::max(lambda / LAMBDA_DN, LAMBDA_MIN); // Diminution du lambda
			}
			else {
				// ACCEPTE quand même (line-search style) car Err_Try < err, mais on reste prudent sur lambda
				err_old = err;
				ws.p[1] = ws.p_try;
				err     = lm(ws);
				lambda  = std::max(lambda / std::sqrt(LAMBDA_DN), LAMBDA_MIN); // baisse plus timide
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

	std::cout << "***** Final Point : " << offset << ", " << intensity << ", " << x << ", " << y << ", " << z << " *****\n";
	fill(infos, point, offset, intensity, x, y, z, 0, 0, 0, rss, rss);
}

//-------------------------------------------------------------------------------------------------
void SplineFit::delta3D(LM_Workspace& ws, const std::vector<double>& p) const
{
	// Coordonnées fractionnaires (ROI locale)
	const double Xr = -1.0 * (p[2] - static_cast<double>(roiHalf_) + 0.5), Yr = -1.0 * (p[3] - static_cast<double>(roiHalf_) + 0.5);

	// Parties entières pour l’ancrage du patch
	const int Xi = static_cast<int>(std::floor(Xr)), Yi = static_cast<int>(std::floor(Yr)), Zi = static_cast<int>(std::floor(p[4]));

	// Parties fractionnaires à l’intérieur de la cellule
	const double Xc = Xr - Xi, Yc = Yr - Yi, Zc = p[4] - Zi;

	auto bspline3 = [](const double u, std::array<double, SPLINE_PATCH_SIDE>& w, std::array<double, SPLINE_PATCH_SIDE>& dw) {
		const double U2 = u * u, U3    = U2 * u;
		const double Om = 1.0 - u, Om2 = Om * Om, Om3 = Om2 * Om;

		// Poids (uniform cubic B-spline)
		w[0] = Om3 / 6.0;
		w[1] = (3.0 * U3 - 6.0 * U2 + 4.0) / 6.0;
		w[2] = (-3.0 * U3 + 3.0 * U2 + 3.0 * u + 1.0) / 6.0;
		w[3] = U3 / 6.0;

		// Dérivées
		dw[0] = -0.5 * Om2;
		dw[1] = 1.5 * U2 - 2.0 * u;
		dw[2] = -1.5 * U2 + u + 0.5;
		dw[3] = 0.5 * U2;
	};

	std::array<double, SPLINE_PATCH_SIDE> wx{}, wy{}, wz{}, dwx{}, dwy{}, dwz{};
	bspline3(Xc, wx, dwx);
	bspline3(Yc, wy, dwy);
	bspline3(Zc, wz, dwz);

	std::ranges::fill(ws.f, 0.0);
	for (auto& df : ws.df) { std::ranges::fill(df, 0.0); }

	for (size_t iz = 0, id = 0; iz < SPLINE_PATCH_SIDE; ++iz) {
		for (size_t iy = 0; iy < SPLINE_PATCH_SIDE; ++iy) {
			for (size_t ix = 0; ix < SPLINE_PATCH_SIDE; ++ix, ++id) {
				const double Wyz = wz[iz] * wy[iy];
				ws.f[id]         = Wyz * wx[ix];
				ws.df[0][id]     = Wyz * dwx[ix];             // ∂s/∂Xc
				ws.df[1][id]     = wz[iz] * dwy[iy] * wx[ix]; // ∂s/∂Yc
				ws.df[2][id]     = dwz[iz] * wy[iy] * wx[ix]; // ∂s/∂Zc
			}
		}
	}

	double swx = 0, swy = 0, swz = 0, sf = 0;
	for (const auto W : wx) { swx += W; }
	for (const auto W : wy) { swy += W; }
	for (const auto W : wz) { swz += W; }
	for (const auto F : ws.f) { sf += F; }
	std::cout << "sum(wx) = " << swx << ", sum(wy) = " << swy << ", sum(wz) = " << swz << ", sum(f) = " << sf << "\n";
}

//-------------------------------------------------------------------------------------------------
double SplineFit::lm(LM_Workspace& ws) const
{
	std::ranges::fill(ws.alpha, 0.0);
	std::ranges::fill(ws.beta, 0.0);
	double err = 0.0;

	delta3D(ws, ws.p[1]); // Mise à jour des contributions d'une base de fonctions splines 3D (valeurs et dérivées).
	// On se base sur les décimales donc f et df mis à jour ici ne changent pas pour toute la ROI
	// (ça à l'air d'etre une sorte de "poids" pour chacun des 64 coeffs)

	// Parties entières pour l’ancrage du patch
	const int Xi = static_cast<int>(std::floor(-1.0 * (ws.p[1][2] - static_cast<double>(roiHalf_) + 0.5))) + static_cast<int>(splineOffset_),
	          Yi = static_cast<int>(std::floor(-1.0 * (ws.p[1][3] - static_cast<double>(roiHalf_) + 0.5))) + static_cast<int>(splineOffset_),
	          Zi = static_cast<int>(std::floor(ws.p[1][4]));

	for (size_t y = 0, idx = 0; y < roi_; ++y) {
		const int Py = static_cast<int>(y) + Yi;
		for (size_t x = 0; x < roi_; ++x, ++idx) {
			const int Px   = static_cast<int>(x) + Xi;
			Vec3I pos      = Vec3I::IntToVec(Px, Py, Zi, splineDim_ - 1);
			ws.hat[idx]    = psf(ws, ws.p[1], pos);          // PSF théorique au pixel (x, y) selon les paramètres p
			ws.delta[idx]  = ws.grid[idx].v - ws.hat[idx];   // Résidu (∂) = mesure - prédiction
			ws.deltaW[idx] = ws.weight[idx] * ws.delta[idx]; // Résidu pondéré (W * ∂)
			err += ws.delta[idx] * ws.deltaW[idx];           // mise à jour du Chi² (somme des carrés pondérée W * ∂ * ∂)
			deriv(ws, pos);                                  // dérivées au pixel

			// Accumulations alpha += w * grad * gradᵀ  et  beta += w * grad * r
			for (size_t a = 0; a < nParams_; ++a) {
				ws.beta[a] -= ws.grad[a] * ws.deltaW[idx];
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
		const double C  = coeffs_[Id] * ws.coeffRatio; // <<< normalisation pour ∂s/∂*
		ws.grad[1] += ws.f[i] * C;                     // Contribution à l’intensité.
		ws.grad[2] += ws.df[0][i] * C;                 // Dérivée en x.
		ws.grad[3] += ws.df[1][i] * C;                 // Dérivée en y.
		ws.grad[4] += ws.df[2][i] * C;                 // Dérivée en z.
	}

	const double I = std::max(ws.p[1][1], ws.Imin); // vraiment utile ? je crois que les différents clamps précédents devraient gérer ça
	ws.grad[0]     = 1;                             // Dérivée par rapport à l’offset
	ws.grad[2] *= -1 * I;                           // Ajustement en x.
	ws.grad[3] *= -1 * I;                           // Ajustement en y.
	ws.grad[4] *= I;                                // Ajustement en z.
}

//-------------------------------------------------------------------------------------------------
double SplineFit::psf(const LM_Workspace& ws, const std::vector<double>& p, const Vec3I& pos) const
{
	const size_t Pos = pos.Clamp(splineDim_ - 1).GetIndex(splineDim_);
	double s         = 0.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) { s += ws.f[i] * coeffs_[i * splineSize_ + Pos]; }
	return p[0] + p[1] * s * ws.coeffRatio; // <<< normalisation avec ws.coeffRatio
}

//-------------------------------------------------------------------------------------------------
double SplineFit::getErr(LM_Workspace& ws) const
{
	double res = 0.0;
	delta3D(ws, ws.p_try); // Mise à jour des contributions d'une base de fonctions splines 3D (valeurs et dérivées).

	// Parties entières pour l’ancrage du patch
	const int Xi = static_cast<int>(std::floor(-1.0 * (ws.p_try[2] - static_cast<double>(roiHalf_) + 0.5))) + static_cast<int>(splineOffset_),
	          Yi = static_cast<int>(std::floor(-1.0 * (ws.p_try[3] - static_cast<double>(roiHalf_) + 0.5))) + static_cast<int>(splineOffset_),
	          Zi = static_cast<int>(std::floor(ws.p_try[4]));

	for (size_t y = 0, idx = 0; y < roi_; ++y) {
		const int Py = static_cast<int>(y) + Yi;
		for (size_t x = 0; x < roi_; ++x, ++idx) {
			const int Px    = static_cast<int>(x) + Xi;
			const Vec3I Pos = Vec3I::IntToVec(Px, Py, Zi, splineDim_ - 1);
			const double Dx = ws.grid[idx].v - psf(ws, ws.p_try, Pos); // PSF théorique au pixel (x, y) selon les paramètres p
			res += ws.weight[idx] * Dx * Dx;
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
void SplineFit::computeCoeffRatio(LM_Workspace& ws, const std::vector<double>& p) const
{
	// Sauvegarde f / df (on ne veut PAS polluer ws pour la suite du LM)
	const auto FSave  = ws.f;
	const auto DfSave = ws.df;

	// Met à jour f au centre (0.5,0.5,0.5)
	delta3D(ws, p);

	// Cellule centrale (incluant l’ancrage entier Xi, Yi du patch courant)
	const int Xi = static_cast<int>(std::floor(-1.0 * (p[2] - static_cast<double>(roiHalf_) + 0.5))) + static_cast<int>(splineOffset_);
	const int Yi = static_cast<int>(std::floor(-1.0 * (p[3] - static_cast<double>(roiHalf_) + 0.5))) + static_cast<int>(splineOffset_);
	const int Zi = static_cast<int>(std::floor(p[4]));

	const size_t Pos = Vec3I::IntToVec(Xi, Yi, Zi, splineDim_ - 1).GetIndex(splineDim_);

	// Somme pondérée au centre de cellule
	double s0 = 0.0, sL1 = 0.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) {
		const double C = coeffs_[i * splineSize_ + Pos],
		             W = ws.f[i]; // ≥ 0, et Σw = 1
		s0 += ws.f[i] * C;
		sL1 += W * std::abs(C);
	}

	// Restaure f / df
	ws.f  = FSave;
	ws.df = DfSave;

	// Garde-fou : impose un |s0| minimal relatif
	const double MinAbs = std::max(1e-3 * sL1, 1e-12);

	// >>> toujours POSITIF
	s0 = std::max(std::abs(s0), MinAbs);

	// Ratio borné, POSITIF
	double ratio          = 1.0 / s0;
	constexpr double Rmax = 1e2; // 1e2 (plus doux que 1e3)
	if (ratio > Rmax) { ratio = Rmax; }

	ws.coeffRatio = ratio;
}
} // namespace palm_tracer
