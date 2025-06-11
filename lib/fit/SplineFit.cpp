#include "SplineFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <iostream>
#include <limits>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "Utils.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void SplineFit::parseParams()
{
	roi_           = static_cast<size_t>(fitParams_[0]);
	splineSizeX_   = static_cast<size_t>(fitParams_[1]);
	splineSizeY_   = static_cast<size_t>(fitParams_[2]);
	splineSizeZ_   = static_cast<size_t>(fitParams_[3]);
	splineSizeXY_  = splineSizeX_ * splineSizeY_;
	splineSizeXYZ_ = splineSizeXY_ * splineSizeZ_;
	initZ_         = fitParams_[4];
	coeffs_.resize(splineSizeXYZ_ * SPLINE_PATCH_SIZE);
	std::ranges::fill(coeffs_, 0.0);
	//for (size_t i = 0; i < coeffs_.size(); ++i) { coeffs_[i] = fitParams_[i + 4]; }
	const double CenterCoeff = coeffs_[splineSizeZ_ / 2 * splineSizeXY_ + splineSizeY_ / 2 * splineSizeX_ + splineSizeX_ / 2];
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
	const size_t Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux <= width_ - roiHalf_ - 1 && roiHalf_ <= Uy && Uy <= height_ - roiHalf_ - 1) {
		// Création des variables
		double lambda = LAMBDA_0, minI = 0, maxI = 0, model = 0, err = std::numeric_limits<double>::max(), mu = 1 + lambda;
		std::vector<double> data(n_), hessian(nParams2_),
		                    dudt(nParams_), jacobian(nParams_), jump(nParams_), update(nParams_), l(nParams_), u(nParams_),
		                    f(SPLINE_PATCH_SIZE), dxf(SPLINE_PATCH_SIZE), dyf(SPLINE_PATCH_SIZE), dzf(SPLINE_PATCH_SIZE);
		std::vector<size_t> coords(2 * n_);
		std::array<std::vector<double>, 3> p;

		// Initialization of data & t (A vérifier)
		preCompute(Ux, Uy, data, coords, minI, maxI);

		// Initialisation du/des P :
		// - Calcul du centre de masse de la région pour x et y
		CenterOfMass(data, roi_, x, y);

		initP(p, x, y, z, minI, maxI);
		p[1][1] *= intensityRatio_; // Change l'intensité pour la position actuelle.
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
			const double Xc = p[1][2] - std::floor(p[1][2]), Yc = p[1][3] - std::floor(p[1][3]), Zc = p[1][4] - std::floor(p[1][4]);
			Delta3D(Xc, Yc, Zc, f, dxf, dyf, dzf);

			// Calcul des dérivées et initialisation du model
			// offset utile, car la ROI et les coeffs ont peut-être des tailles différentes donc on recentre sur X et sur Y
			// Exemple : ROI de 7 splineX de 11 : La ROI 7x7 est centrée à partir du point 2 dans la spline
			for (size_t i = 0; i < n_; ++i) {
				derivativeSpline(coords[i * 2] + splineOffset_, coords[i * 2 + 1] + splineOffset_, static_cast<size_t>(std::floor(p[1][4])),
				                 f, dxf, dyf, dzf, p[1], dudt, model);
				updateModel(data[i], model, err, jacobian, hessian, dudt);
			}

			// Mise à jour de la matrice Hessian
			for (size_t i = 0; i < nParams_; ++i) { hessian[i * nParams_ + i] = hessian[i * nParams_ + i] * mu; }

			//----- Changements de P, mu, lambda -----
			if (Cholesky(hessian, nParams_, l, u)) { // Si la Factorisation de Cholesky a réussi
				// Calcul du vecteur de mise à jour
				luEvaluate(l, u, jacobian, update);
				// Mise à jour de la matrice de saut, la matrice update et de P.
				for (size_t i = 0; i < nParams_; ++i) {
					if (update[i] / OldUpdate[i] < -0.5) { jump[i] *= 0.5; }
					update[i] /= (1 + fabs(update[i] / jump[i]));
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

		offset    = p[1][0];
		intensity = p[1][1];
		x         = p[1][2];
		y         = p[1][3];
		z         = p[1][4];
		//Calculs du MSE sur XY et sur Z à faire...
	}

	fill(infos, point, offset, intensity, x, y, z, 0, 0, 0, rssXY, rssZ);
}

//-------------------------------------------------------------------------------------------------
void SplineFit::derivativeSpline(const size_t x, const size_t y, const size_t z, const std::vector<double>& f,
                                 const std::vector<double>& dxf, const std::vector<double>& dyf, const std::vector<double>& dzf,
                                 const std::vector<double>& p, std::vector<double>& dudt, double& model) const
{
	const size_t Pos = z * splineSizeXY_ + y * splineSizeX_ + x;

	for (size_t i = 0; i < SPLINE_PATCH_SIZE; ++i) {
		const size_t Id = i * splineSizeXYZ_ + Pos;
		dudt[1] += f[i] * coeffs_[Id];   // Intensité
		dudt[2] += dxf[i] * coeffs_[Id]; // x
		dudt[3] += dyf[i] * coeffs_[Id]; // y
		dudt[4] += dzf[i] * coeffs_[Id]; // z
	}

	dudt[0] = 1;          // offset
	dudt[2] *= -1 * p[1]; // x
	dudt[3] *= -1 * p[1]; // y
	dudt[4] *= p[1];      // z
	model = p[0] + p[1] * dudt[1];
}

void SplineFit::updateModel(const double data, const double model, double& err, std::vector<double>& jacobian,
                            std::vector<double>& hessian, const std::vector<double>& dudt) const
{
	if (data > 0) { err = err + 2 * (model - data - data * log(model / data)); }
	else { err = err + 2 * model; }

	const double T1 = 1 - data / model;
	for (size_t i = 0; i < nParams_; ++i) { jacobian[i] += T1 * dudt[i]; }

	const double T2 = data / (model * model);
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = i; j < nParams_; ++j) {
			hessian[i * nParams_ + j] += T2 * dudt[i] * dudt[j];
			hessian[j * nParams_ + i] = hessian[i * nParams_ + j];
		}
	}
}

void SplineFit::updateMuLambda(double& mu, double& lambda, const bool up)
{
	if (up) {
		mu     = std::max((1 + lambda * LAMBDA_UP) / (1 + lambda), 1.3);
		lambda = LAMBDA_UP * lambda;
	}
	else {
		lambda = LAMBDA_DN * lambda;
		mu     = 1 + lambda;
	}
}

void SplineFit::luEvaluate(const std::vector<double>& l, const std::vector<double>& u, const std::vector<double>& b, std::vector<double>& x) const
{
	//ax = b -> lux = b. Then y is defined to be Ux
	std::vector y(6, 0.0);
	// Forward solve Ly = b
	for (size_t i = 0; i < nParams_; ++i) {
		y[i] = b[i];
		for (size_t j = 0; j < i; ++j) { y[i] -= l[j * nParams_ + i] * y[j]; }
		y[i] /= l[i * nParams_ + i];
	}
	// Backward solve Ux = y
	for (int i = static_cast<int>(nParams_ - 1); i >= 0; --i) {
		x[i] = y[i];
		for (size_t j = i + 1; j < nParams_; j++) { x[i] -= u[j * nParams_ + i] * x[j]; }
		x[i] /= u[i * nParams_ + i];
	}
}
}
