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

//==================== Common Functions ====================
//-------------------------------------------------------------------------------------------------
void SplineFit::parseParams()
{
	size_t idx = 0;
	roi_       = static_cast<size_t>(fitParams_[idx++]);
	for (size_t i = 0; i < 3; ++i) { splineDim_[i] = static_cast<size_t>(fitParams_[idx++]); }
	splineSizeXY_ = splineDim_.GetPlaneSize();
	splineSize_   = splineDim_.GetSize();
	initZ_        = static_cast<double>(splineDim_.z) / 2.0;
	deltaZ_       = fitParams_[idx++];
	coeffs_.resize(splineSize_ * SPLINE_PATCH_SIZE);
	for (double& coeff : coeffs_) { coeff = fitParams_[idx++]; }
	const double CenterCoeff = coeffs_[32 * splineSize_ + splineDim_.z / 2 * splineSizeXY_ + splineDim_.y / 2 * splineDim_.x + splineDim_.x / 2];
	intensityRatio_          = std::abs(CenterCoeff) < std::numeric_limits<double>::epsilon() ? 1.0 : 4.0 / CenterCoeff;
	splineOffset_            = static_cast<size_t>(std::floor((static_cast<double>(splineDim_.x) + 1.0 - static_cast<double>(roi_)) / 2));
}

//-------------------------------------------------------------------------------------------------
void SplineFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const
{
	double offset    = 0,
	       intensity = 0,
	       x         = infos.x_,
	       y         = infos.y_,
	       z         = initZ_,
	       rssXY     = -1,
	       rssZ      = -1;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const auto Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux < dim_.x - roiHalf_ && roiHalf_ <= Uy && Uy < dim_.y - roiHalf_) {
		// Création des variables
		double minI = 0, maxI = 0;
		MLE_Workspace ws;
		ws.Resize(n_, nParams_, nParams2_);

		// Initialization of data & t (A vérifier)
		preCompute(Ux, Uy, ws.grid, minI, maxI);

		// Initialisation du/des P :
		// - Calcul du centre de masse de la région pour x et y
		matrix::CenterOfMass(ws.grid, roi_, x, y);

		initP(ws.p, x, y, z, minI, maxI);
		ws.p[1][1] *= intensityRatio_;                                                                         // Change l'intensité pour la position actuelle.
		ws.p[2][4] = static_cast<double>(splineDim_.z);                                                        // Change le maximum pour Z.
		for (size_t i = 0; i < nParams_; ++i) { ws.p[1][i] = std::clamp(ws.p[1][i], ws.p[0][i], ws.p[2][i]); } // Clamp en cas de problème.
		ws.jump = {std::max(20.0, ws.p[1][1]), std::max(100.0, ws.p[1][1]), 1.0, 1.0, std::max(static_cast<double>(splineDim_.z) / 3.0, 2.0)};

		//--------------------------------------------------------------------------------------------------------------------------------------------
		std::cout << std::setprecision(6) << std::fixed;
		std::cout << "\t [min .. val .. max]\n";
		constexpr std::array PName{"B", "I", "X", "Y", "Z"};
		for (size_t i = 0; i < nParams_; ++i) { std::cout << "\t" << PName[i] << "= [" << ws.p[0][i] << " .. " << ws.p[1][i] << " .. " << ws.p[2][i] << "]\n"; }
		std::cout << "Jump";
		matrix::PrintVector(ws.jump);
		std::cout << "\n";
		//--------------------------------------------------------------------------------------------------------------------------------------------

		bool stop        = false;
		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			//--------------------------------------------------------------------------------------------------------------------------------------------
			std::cout << "\n--------------------------------------------------\n\t Iteration " << iteration << "\n";
			//--------------------------------------------------------------------------------------------------------------------------------------------

			//----- Conservation des anciennes valeurs de l'erreur, du point, de l'update, et du lambda -----
			const double OldErr                 = ws.err;
			const std::vector<double> OldP      = ws.p[1];
			const std::vector<double> OldUpdate = ws.update;
			const double OldLambda              = ws.lambda;
			//----- Remise à zéro des vecteurs -----
			std::ranges::fill(ws.jacobian, 0.0);
			std::ranges::fill(ws.hessian, 0.0);
			std::ranges::fill(ws.update, 0.0);
			std::ranges::fill(ws.l, 0.0);
			std::ranges::fill(ws.u, 0.0);

			//----- Calculs -----
			// Calcul du patch 3D de 4x4x4 (64) (SPLINE_PATCH_SIDE = 4, SPLINE_PATCH_SIZE = 64)
			const Vec3D PosR{
				.x = -1 * (ws.p[1][2] - static_cast<double>(roiHalf_) + 0.5),
				.y = -1 * (ws.p[1][3] - static_cast<double>(roiHalf_) + 0.5),
				.z = ws.p[1][4]
			};

			const int Xi = static_cast<int>(std::floor(PosR.x)),
			          Yi = static_cast<int>(std::floor(PosR.y)),
			          Zi = static_cast<int>(std::floor(PosR.z));

			const Vec3D PosC{.x = PosR.x - Xi, .y = PosR.y - Yi, .z = PosR.z - Zi};
			delta3D(PosC, ws);
			// Calcul des dérivées et initialisation du model
			// offset utile, car la ROI et les coeffs ont peut-être des tailles différentes donc on recentre sur X et sur Y
			// Exemple : ROI de 7 splineX de 11 : La ROI 7x7 est centrée à partir du point 2 dans la spline
			Vec3I pos = {.x = 0, .y = 0, .z = static_cast<size_t>(Zi)};
			for (size_t i = 0, idx = 0; i < roi_; ++i) {
				pos.x = Xi + i + splineOffset_;
				for (size_t j = 0; j < roi_; ++j, ++idx) {
					pos.y = Yi + j + splineOffset_;
					derivativeSpline(pos, ws);
					updateJacobianHessian(ws.grid[idx].v, ws);
				}
			}

			// Mise à jour de la matrice Hessian
			for (size_t i = 0; i < nParams_; ++i) { ws.hessian[i * nParams_ + i] = ws.hessian[i * nParams_ + i] * ws.mu; }

			//--------------------------------------------------------------------------------------------------------------------------------------------
			std::cout << "\tJacobian :";
			matrix::PrintVector(ws.jacobian);
			std::cout << "\n";
			std::cout << "\tHessian :\n";
			matrix::PrintMatrix(ws.hessian, nParams_, nParams_);
			std::cout << "\n";

			std::cout << "\tPoint (" << ws.p[1][2] << ", " << ws.p[1][3] << ", " << ws.p[1][4]
				<< ")\n\terr = " << ws.err << "\t lambda = " << ws.lambda << "\t model = " << ws.model << "\n";
			//--------------------------------------------------------------------------------------------------------------------------------------------

			//----- Changements de P, mu, lambda -----
			if (matrix::Cholesky(ws.hessian, nParams_, ws.l, ws.u)) { // Si la Factorisation de Cholesky a réussi
				// Calcul du vecteur de mise à jour
				luEvaluate(ws);
				// Mise à jour de la matrice de saut, la matrice update et de P.
				for (size_t i = 0; i < nParams_; ++i) {
					if (std::abs(OldUpdate[i]) > EPSILON && ws.update[i] / OldUpdate[i] < -0.5) { ws.jump[i] *= 0.5; }
					const double Scale = ws.update[i] / (1.0 + std::abs(ws.update[i]) / ws.jump[i]); // Pas de division par une grandeur instable
					ws.p[1][i] -= Scale;
					ws.p[1][i]   = std::clamp(ws.p[1][i], ws.p[0][i], ws.p[2][i]);
					ws.update[i] = Scale;
				}
				if (ws.err > OldErr) { // si erreur plus importante, restauration des valeurs et augmentation du lambda
					ws.err    = OldErr;
					ws.p[1]   = OldP;
					ws.update = OldUpdate;
					ws.lambda = OldLambda;
					updateMuLambda(ws, true);
				}
				else { updateMuLambda(ws, false); } // sinon réduction
				std::cout << "\t\tAfter Update lambda = " << ws.lambda << "\t mu = " << ws.mu << "\n";
				std::cout << "\t\tL :\n";
				matrix::PrintMatrix(ws.l, nParams_, nParams_);
				std::cout << "\n";
				std::cout << "\t\tU :\n";
				matrix::PrintMatrix(ws.l, nParams_, nParams_);
				std::cout << "\n";
				std::cout << "\t\tUpdate :";
				matrix::PrintVector(ws.update);
				std::cout << "\n";
			}
			else { // restauration des valeurs et Augmentation du lambda si la factorisation a échoué
				ws.err    = OldErr;
				ws.p[1]   = OldP;
				ws.update = OldUpdate;
				ws.lambda = OldLambda;
				updateMuLambda(ws, true);
			}

			//----- Conditions d'arrêts -----
			if (iteration >= N_ITER) { stop = true; } // Si Nombre d'itérations maximum atteint, on arrête.
			if (iteration > 2 && (ws.err < EPSILON_1 || std::fabs(ws.err - OldErr) / ws.err < EPSILON_1)) { stop = true; } // Si Convergence, on arrête.
		}

		//Calcul du CRLB
		std::vector<double> m(nParams2_, 0.0), mInv(nParams2_, 0.0), crlb(nParams_);
		std::ranges::fill(ws.dudt, 0.0);

		const Vec3D PosR{
			.x = -1 * (ws.p[1][2] - static_cast<double>(roiHalf_) + 0.5),
			.y = -1 * (ws.p[1][3] - static_cast<double>(roiHalf_) + 0.5),
			.z = ws.p[1][4]
		};

		const int Xi = static_cast<int>(std::floor(PosR.x)),
		          Yi = static_cast<int>(std::floor(PosR.y)),
		          Zi = static_cast<int>(std::floor(PosR.z));

		const Vec3D PosC{.x = PosR.x - Xi, .y = PosR.y - Yi, .z = PosR.z - Zi};
		delta3D(PosC, ws);
		// Calcul des dérivées et initialisation du model
		// offset utile, car la ROI et les coeffs ont peut-être des tailles différentes donc on recentre sur X et sur Y
		// Exemple : ROI de 7 splineX de 11 : La ROI 7x7 est centrée à partir du point 2 dans la spline
		Vec3I pos = {.x = 0, .y = 0, .z = static_cast<size_t>(Zi)};
		for (size_t i = 0, idx = 0; i < roi_; ++i) {
			pos.x = Xi + i + splineOffset_;
			for (size_t j = 0; j < roi_; ++j, ++idx) {
				pos.y = Yi + j + splineOffset_;
				derivativeSpline(pos, ws);
				// Calcul de la matrice d'information de Fisher
				for (size_t i2 = 0; i2 < nParams_; ++i2) {
					for (size_t j2 = i2; j2 < nParams_; ++j2) {
						m[i2 * nParams_ + j2] += ws.dudt[j2] * ws.dudt[i2] / ws.model;
						m[j2 * nParams_ + i2] = m[i2 * nParams_ + j2];
					}
				}
			}
		}


		// Si la convergence est parfaite (ou le point pile sur une grille), des 0 seront présents, à remplacer pour calculer la matrice inverse.
		for (size_t i = 0; i < nParams_; ++i) {
			if (std::abs(m[i * nParams_ + i]) < std::numeric_limits<double>::epsilon()) { m[i * nParams_ + i] = EPSILON_2; }
		}

		// Matrice inverse (CRLB=F^-1)
		matrix::InverseLU(m, mInv, nParams_);

		for (size_t i = 0; i < nParams_; ++i) { crlb[i] = SqrtReal(mInv[i * nParams_ + i]); }

		// Assignation
		offset    = crlb[0];
		intensity = crlb[1];
		x         = ws.p[1][2] - static_cast<double>(roiHalf_) + infos.x_;
		y         = ws.p[1][3] - static_cast<double>(roiHalf_) + infos.y_;
		z         = deltaZ_ * (ws.p[1][4] - initZ_); //ws.p[1][4]; //deltaZ_ * (ws.p[1][4] - initZ_);
		rssXY     = (crlb[2] + crlb[3]) / 2;
		rssZ      = deltaZ_ * crlb[4];
	}
	else { z = 0; }

	fill(infos, point, offset, intensity, x, y, z, 0, 0, 0, rssXY, rssZ);
}

//==================== Maximum de vraisemblance ====================
//-------------------------------------------------------------------------------------------------
void SplineFit::updateJacobianHessian(const double data, MLE_Workspace& ws) const
{
	// Calcul du terme d’erreur : -2 * log-vraisemblance (formule EMCCD)
	if (data > 0) { ws.err = ws.err + 2 * (ws.model - data - data * log(ws.model / data)); }
	else { ws.err = ws.err + 2 * ws.model; }

	const double T1 = 1 - data / ws.model;
	for (size_t i = 0; i < nParams_; ++i) { ws.jacobian[i] += T1 * ws.dudt[i]; } // ∂L/∂θ_i

	const double T2 = data / (ws.model * ws.model);
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = i; j < nParams_; ++j) {
			ws.hessian[i * nParams_ + j] += T2 * ws.dudt[i] * ws.dudt[j];
			ws.hessian[j * nParams_ + i] = ws.hessian[i * nParams_ + j]; // Hessienne symétrique
		}
	}
}

//-------------------------------------------------------------------------------------------------
void SplineFit::updateMuLambda(MLE_Workspace& ws, const bool up)
{
	if (up) { // Augmentation du lambda : on recule dans la descente
		ws.mu = std::max((1 + ws.lambda * LAMBDA_UP) / (1 + ws.lambda), 1.3);
		ws.lambda *= LAMBDA_UP;
	}
	else { // Diminution du lambda : on fait confiance au gradient
		ws.lambda /= LAMBDA_DN;
		ws.mu = 1 + ws.lambda;
	}
}

//-------------------------------------------------------------------------------------------------
void SplineFit::luEvaluate(MLE_Workspace& ws) const
{
	// ax = b ⇾ lux = b, on résout Ly = b puis Ux = y
	std::vector y(nParams_, 0.0);

	// Substitution avant : Ly = b
	for (size_t i = 0; i < nParams_; ++i) {
		y[i] = ws.jacobian[i];
		for (size_t j = 0; j < i; ++j) { y[i] -= ws.l[j * nParams_ + i] * y[j]; }
		y[i] /= ws.l[i * nParams_ + i];
	}

	// Substitution arrière : Ux = y
	for (int i = static_cast<int>(nParams_ - 1); i >= 0; --i) {
		ws.update[i] = y[i];
		for (size_t j = i + 1; j < nParams_; j++) { ws.update[i] -= ws.u[j * nParams_ + i] * ws.update[j]; }
		ws.update[i] /= ws.u[i * nParams_ + i];
	}
}

//==================== Run Computes ====================
//-------------------------------------------------------------------------------------------------
void SplineFit::derivativeSpline(Vec3I& pos, MLE_Workspace& ws) const
{
	pos = pos.Clamp(splineDim_ - 1);

	const size_t Pos = pos.GetIndex(splineDim_);
	std::ranges::fill(ws.dudt, 0.0); // Initialise toutes les dérivées à zéro


	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) {
		const size_t Id = i * splineSize_ + Pos;
		ws.dudt[1] += ws.f[i] * coeffs_[Id];     // Contribution à l’intensité.
		ws.dudt[2] += ws.df[0][i] * coeffs_[Id]; // Dérivée en x.
		ws.dudt[3] += ws.df[1][i] * coeffs_[Id]; // Dérivée en y.
		ws.dudt[4] += ws.df[2][i] * coeffs_[Id]; // Dérivée en z.
	}
	ws.dudt[0] = 1;                                  // Dérivée par rapport à l’offset
	ws.dudt[2] *= -1 * ws.p[1][1];                   // Ajustement en x.
	ws.dudt[3] *= -1 * ws.p[1][1];                   // Ajustement en y.
	ws.dudt[4] *= ws.p[1][1];                        // Ajustement en z.
	ws.model = ws.p[1][0] + ws.p[1][1] * ws.dudt[1]; // Évaluation du modèle spline : offset + amplitude * contribution
}

//-------------------------------------------------------------------------------------------------
void SplineFit::delta3D(const Vec3D& pos, MLE_Workspace& ws)
{
	std::ranges::fill(ws.f, 0.0);
	for (size_t i = 0; i < 3; ++i) { std::ranges::fill(ws.df[i], 0.0); }
	double cz = 1.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIDE; ++i) {
		double cy = 1.0;
		for (size_t j = 0; j < SPLINE_PATCH_SIDE; ++j) {
			double cx = 1.0;
			for (size_t k = 0; k < SPLINE_PATCH_SIDE; ++k) {
				const size_t Id   = i * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k;
				const double Coef = cz * cy * cx;
				ws.f[Id]          = Coef;
				if (k + 1 < SPLINE_PATCH_SIDE) { ws.df[0][Id + 1] = static_cast<double>(k + 1) * Coef; }
				if (j + 1 < SPLINE_PATCH_SIDE) { ws.df[1][Id + SPLINE_PATCH_SIDE] = static_cast<double>(j + 1) * Coef; }
				if (i + 1 < SPLINE_PATCH_SIDE) { ws.df[2][Id + SPLINE_PATCH_PLANE] = static_cast<double>(i + 1) * Coef; }
				cx *= pos.x;
			}
			cy *= pos.y;
		}
		cz *= pos.z;
	}
}
} // namespace palm_tracer
