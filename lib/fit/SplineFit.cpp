#include "SplineFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <execution>
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
	size_t idx     = 0;
	roi_           = static_cast<size_t>(fitParams_[idx++]);
	splineSizeX_   = static_cast<size_t>(fitParams_[idx++]);
	splineSizeY_   = static_cast<size_t>(fitParams_[idx++]);
	splineSizeZ_   = static_cast<size_t>(fitParams_[idx++]);
	splineSizeXY_  = splineSizeX_ * splineSizeY_;
	splineSizeXYZ_ = splineSizeXY_ * splineSizeZ_;
	initZ_         = static_cast<double>(splineSizeZ_) / 2.0;
	deltaZ_        = fitParams_[idx++];
	coeffs_.resize(splineSizeXYZ_ * SPLINE_PATCH_SIZE);
	for (double& coeff : coeffs_) { coeff = fitParams_[idx++]; }
	const double CenterCoeff = coeffs_[32 * splineSizeXYZ_ + splineSizeZ_ / 2 * splineSizeXY_ + splineSizeY_ / 2 * splineSizeX_ + splineSizeX_ / 2];
	intensityRatio_          = std::fabs(CenterCoeff) < std::numeric_limits<double>::epsilon() ? 1.0 : 4.0 / CenterCoeff;
	splineOffset_            = static_cast<size_t>(std::floor((static_cast<double>(splineSizeX_) + 1.0 - static_cast<double>(roi_)) / 2));
}

//-------------------------------------------------------------------------------------------------
void SplineFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point)
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
	if (roiHalf_ <= Ux && Ux < width_ - roiHalf_ && roiHalf_ <= Uy && Uy < height_ - roiHalf_) {
		// Création des variables
		double lambda = LAMBDA_0, minI = 0, maxI = 0, model = 0, err = std::numeric_limits<double>::max(), mu = 1 + lambda;
		std::vector<double> hessian(nParams2_), l(nParams2_), u(nParams2_),
		                    dudt(nParams_), jacobian(nParams_), jump(nParams_), update(nParams_),
		                    f(SPLINE_PATCH_SIZE), dxf(SPLINE_PATCH_SIZE), dyf(SPLINE_PATCH_SIZE), dzf(SPLINE_PATCH_SIZE);
		std::vector<Pixel> grid(n_);
		std::array<std::vector<double>, 3> p;

		// Initialization of data & t (A vérifier)
		preCompute(Ux, Uy, grid, minI, maxI);

		// Initialisation du/des P :
		// - Calcul du centre de masse de la région pour x et y
		matrix::CenterOfMass(grid, roi_, x, y);

		initP(p, x, y, z, minI, maxI);
		p[1][1] *= intensityRatio_;                                                                // Change l'intensité pour la position actuelle.
		p[2][4] = static_cast<double>(splineSizeZ_);                                               // Change le maximum pour Z.
		for (size_t i = 0; i < nParams_; ++i) { p[1][i] = std::clamp(p[1][i], p[0][i], p[2][i]); } // Clamp en cas de problème.
		jump = {std::max(20.0, p[1][1]), std::max(100.0, p[1][1]), 1.0, 1.0, std::max(static_cast<double>(splineSizeZ_) / 3.0, 2.0)};

		bool stop        = false;
		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			//----- Conservation des anciennes valeurs de l'erreur, du point, de l'update, et du lambda -----
			const double OldErr                 = err;
			const std::vector<double> OldP      = p[1];
			const std::vector<double> OldUpdate = update;
			const double OldLambda              = lambda;
			//----- Remise à zéro des vecteurs -----
			std::ranges::fill(f, 0.0);
			std::ranges::fill(dxf, 0.0);
			std::ranges::fill(dyf, 0.0);
			std::ranges::fill(dzf, 0.0);
			std::ranges::fill(dudt, 0.0);
			std::ranges::fill(jacobian, 0.0);
			std::ranges::fill(hessian, 0.0);
			std::ranges::fill(update, 0.0);
			std::ranges::fill(l, 0.0);
			std::ranges::fill(u, 0.0);

			//----- Calculs -----
			// Calcul du patch 3D de 4x4x4 (64) (SPLINE_PATCH_SIDE = 4, SPLINE_PATCH_SIZE = 64)
			const double Xr = -1 * (p[1][2] - static_cast<double>(roiHalf_) + 0.5),
			             Yr = -1 * (p[1][3] - static_cast<double>(roiHalf_) + 0.5);

			const int Xi = static_cast<int>(std::floor(Xr)),
			          Yi = static_cast<int>(std::floor(Yr)),
			          Zi = static_cast<int>(std::floor(p[1][4]));

			const double Xc = Xr - Xi, Yc = Yr - Yi, Zc = p[1][4] - Zi;
			delta3D(Xc, Yc, Zc, f, dxf, dyf, dzf);
			// Calcul des dérivées et initialisation du model
			// offset utile, car la ROI et les coeffs ont peut-être des tailles différentes donc on recentre sur X et sur Y
			// Exemple : ROI de 7 splineX de 11 : La ROI 7x7 est centrée à partir du point 2 dans la spline
			for (size_t i = 0, idx = 0; i < roi_; ++i) {
				for (size_t j = 0; j < roi_; ++j, ++idx) {
					derivativeSpline(Xi + i + splineOffset_, Yi + j + splineOffset_, Zi, f, dxf, dyf, dzf, p[1], dudt, model);
					updateJacobianHessian(grid[idx].v, model, err, jacobian, hessian, dudt);
				}
			}

			// Mise à jour de la matrice Hessian
			for (size_t i = 0; i < nParams_; ++i) { hessian[i * nParams_ + i] = hessian[i * nParams_ + i] * mu; }

			//----- Changements de P, mu, lambda -----
			if (matrix::Cholesky(hessian, nParams_, l, u)) { // Si la Factorisation de Cholesky a réussi
				// Calcul du vecteur de mise à jour
				luEvaluate(l, u, jacobian, update);
				// Mise à jour de la matrice de saut, la matrice update et de P.
				for (size_t i = 0; i < nParams_; ++i) {
					if (std::abs(OldUpdate[i]) > EPSILON) {
						if (update[i] / OldUpdate[i] < -0.5) { jump[i] *= 0.5; }
						update[i] /= 1 + std::abs(update[i] / jump[i]);
					}
					else { update[i] /= 1 + std::abs(update[i]) / jump[i]; } // Pas de division par une grandeur instable
					p[1][i] -= update[i];
					p[1][i] = std::clamp(p[1][i], p[0][i], p[2][i]);
				}
				if (err > OldErr) { // si erreur plus importante, restauration des valeurs et augmentation du lambda
					err    = OldErr;
					p[1]   = OldP;
					update = OldUpdate;
					lambda = OldLambda;
					updateMuLambda(mu, lambda, true);
				}
				else { updateMuLambda(mu, lambda, false); } // sinon réduction
			}
			else { // restauration des valeurs et Augmentation du lambda si la factorisation a échoué
				err    = OldErr;
				p[1]   = OldP;
				update = OldUpdate;
				lambda = OldLambda;
				updateMuLambda(mu, lambda, true);
			}

			//----- Conditions d'arrêts -----
			if (iteration >= N_ITER) { stop = true; } // Si Nombre d'itérations maximum atteint, on arrête.
			if (iteration > 2 && (err < EPSILON_1 || std::fabs(err - OldErr) / err < EPSILON_1)) { stop = true; } // Si Convergence, on arrête.
		}

		//Calcul du CRLB
		std::vector<double> m(nParams2_, 0.0), mInv(nParams2_, 0.0), crlb(nParams_);
		std::ranges::fill(f, 0.0);
		std::ranges::fill(dxf, 0.0);
		std::ranges::fill(dyf, 0.0);
		std::ranges::fill(dzf, 0.0);
		std::ranges::fill(dudt, 0.0);

		const double Xr = -1 * (p[1][2] - static_cast<double>(roiHalf_) + 0.5),
		             Yr = -1 * (p[1][3] - static_cast<double>(roiHalf_) + 0.5);

		const int Xi = static_cast<int>(std::floor(Xr)),
		          Yi = static_cast<int>(std::floor(Yr)),
		          Zi = static_cast<int>(std::floor(p[1][4]));

		const double Xc = Xr - Xi, Yc = Yr - Yi, Zc = p[1][4] - Zi;
		delta3D(Xc, Yc, Zc, f, dxf, dyf, dzf);
		// Calcul des dérivées et initialisation du model
		// offset utile, car la ROI et les coeffs ont peut-être des tailles différentes donc on recentre sur X et sur Y
		// Exemple : ROI de 7 splineX de 11 : La ROI 7x7 est centrée à partir du point 2 dans la spline
		for (size_t i = 0, idx = 0; i < roi_; ++i) {
			for (size_t j = 0; j < roi_; ++j, ++idx) {
				derivativeSpline(Xi + i + splineOffset_, Yi + j + splineOffset_, Zi, f, dxf, dyf, dzf, p[1], dudt, model);
				// Calcul de la matrice d'information de Fisher
				for (size_t i2 = 0; i2 < nParams_; ++i2) {
					for (size_t j2 = i2; j2 < nParams_; ++j2) {
						m[i2 * nParams_ + j2] += dudt[j2] * dudt[i2] / model;
						m[j2 * nParams_ + i2] = m[i2 * nParams_ + j2];
					}
				}
			}
		}


		// Si la convergence est parfaite (ou le point pile sur une grille), des 0 seront présents, à remplacer pour calculer la matrice inverse.
		for (size_t i = 0; i < nParams_; ++i) {
			if (std::fabs(m[i * nParams_ + i]) < std::numeric_limits<double>::epsilon()) { m[i * nParams_ + i] = EPSILON_2; }
		}

		// Matrice inverse (CRLB=F^-1)
		matrix::InverseLU(m, mInv, nParams_);

		for (size_t i = 0; i < nParams_; ++i) { crlb[i] = SqrtReal(mInv[i * nParams_ + i]); }

		// Assignation
		offset    = crlb[0];
		intensity = crlb[1];
		x         = p[1][2] - static_cast<double>(roiHalf_) + infos.x_;
		y         = p[1][3] - static_cast<double>(roiHalf_) + infos.y_;
		z         = deltaZ_ * (p[1][4] - initZ_); //p[1][4]; //deltaZ_ * (p[1][4] - initZ_);
		rssXY     = (crlb[2] + crlb[3]) / 2;
		rssZ      = deltaZ_ * crlb[4];
	}
	else { z = 0; }

	fill(infos, point, offset, intensity, x, y, z, 0, 0, 0, rssXY, rssZ);
}

//-------------------------------------------------------------------------------------------------
void SplineFit::derivativeSpline(size_t x, size_t y, size_t z, const std::vector<double>& f,
                                 const std::vector<double>& dxf, const std::vector<double>& dyf, const std::vector<double>& dzf,
                                 const std::vector<double>& p, std::vector<double>& dudt, double& model) const
{
	x = std::clamp<size_t>(x, 0, splineSizeX_ - 1);
	y = std::clamp<size_t>(y, 0, splineSizeY_ - 1);
	z = std::clamp<size_t>(z, 0, splineSizeZ_ - 1);

	const size_t Pos = z * splineSizeXY_ + y * splineSizeX_ + x;
	std::ranges::fill(dudt, 0.0); // Initialise toutes les dérivées à zéro


	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) {
		const size_t Id = i * splineSizeXYZ_ + Pos;
		dudt[1] += f[i] * coeffs_[Id];   // Contribution à l’intensité.
		dudt[2] += dxf[i] * coeffs_[Id]; // Dérivée en x.
		dudt[3] += dyf[i] * coeffs_[Id]; // Dérivée en y.
		dudt[4] += dzf[i] * coeffs_[Id]; // Dérivée en z.
	}
	dudt[0] = 1;                   // Dérivée par rapport à l’offset
	dudt[2] *= -1 * p[1];          // Ajustement en x.
	dudt[3] *= -1 * p[1];          // Ajustement en y.
	dudt[4] *= p[1];               // Ajustement en z.
	model = p[0] + p[1] * dudt[1]; // Évaluation du modèle spline : offset + amplitude * contribution
}

//-------------------------------------------------------------------------------------------------
void SplineFit::updateJacobianHessian(const double data, const double model, double& err, std::vector<double>& jacobian,
                                      std::vector<double>& hessian, const std::vector<double>& dudt) const
{
	// Calcul du terme d’erreur : -2 * log-vraisemblance (formule EMCCD)
	if (data > 0) { err = err + 2 * (model - data - data * log(model / data)); }
	else { err = err + 2 * model; }

	const double T1 = 1 - data / model;
	for (size_t i = 0; i < nParams_; ++i) { jacobian[i] += T1 * dudt[i]; } // ∂L/∂θ_i

	const double T2 = data / (model * model);
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = i; j < nParams_; ++j) {
			hessian[i * nParams_ + j] += T2 * dudt[i] * dudt[j];
			hessian[j * nParams_ + i] = hessian[i * nParams_ + j]; // Hessienne symétrique
		}
	}
}

//-------------------------------------------------------------------------------------------------
void SplineFit::updateMuLambda(double& mu, double& lambda, const bool up)
{
	if (up) { // Augmentation du lambda : on recule dans la descente
		mu     = std::max((1 + lambda * LAMBDA_UP) / (1 + lambda), 1.3);
		lambda = LAMBDA_UP * lambda;
	}
	else { // Diminution du lambda : on fait confiance au gradient
		lambda = LAMBDA_DN * lambda;
		mu     = 1 + lambda;
	}
}

//-------------------------------------------------------------------------------------------------
void SplineFit::luEvaluate(const std::vector<double>& l, const std::vector<double>& u, const std::vector<double>& b, std::vector<double>& x) const
{
	// ax = b ⇾ lux = b, on résout Ly = b puis Ux = y
	std::vector y(nParams_, 0.0);

	// Substitution avant : Ly = b
	for (size_t i = 0; i < nParams_; ++i) {
		y[i] = b[i];
		for (size_t j = 0; j < i; ++j) { y[i] -= l[j * nParams_ + i] * y[j]; }
		y[i] /= l[i * nParams_ + i];
	}

	// Substitution arrière : Ux = y
	for (int i = static_cast<int>(nParams_ - 1); i >= 0; --i) {
		x[i] = y[i];
		for (size_t j = i + 1; j < nParams_; j++) { x[i] -= u[j * nParams_ + i] * x[j]; }
		x[i] /= u[i * nParams_ + i];
	}
}

void SplineFit::delta3D(const double x, const double y, const double z, std::vector<double>& f, std::vector<double>& dxf, std::vector<double>& dyf,
                        std::vector<double>& dzf)
{
	double cz = 1.0;
	for (size_t i = 0; i < SPLINE_PATCH_SIDE; ++i) {
		double cy = 1.0;
		for (size_t j = 0; j < SPLINE_PATCH_SIDE; ++j) {
			double cx = 1.0;
			for (size_t k = 0; k < SPLINE_PATCH_SIDE; ++k) {
				const double Coef                                     = cz * cy * cx;
				f[i * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k] = Coef;
				if (k < SPLINE_PATCH_SIDE - 1) { dxf[i * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k + 1] = static_cast<double>(k + 1) * Coef; }
				if (j < SPLINE_PATCH_SIDE - 1) { dyf[i * SPLINE_PATCH_PLANE + (j + 1) * SPLINE_PATCH_SIDE + k] = static_cast<double>(j + 1) * Coef; }
				if (i < SPLINE_PATCH_SIDE - 1) { dzf[(i + 1) * SPLINE_PATCH_PLANE + j * SPLINE_PATCH_SIDE + k] = static_cast<double>(i + 1) * Coef; }
				cx *= x;
			}
			cy *= y;
		}
		cz *= z;
	}
}
} // namespace palm_tracer
