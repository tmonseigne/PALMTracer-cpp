#include "SplineFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
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
	coeffs_.resize(static_cast<size_t>(fitParams_[5]));
	for (size_t i = 0; i < coeffs_.size(); ++i) { coeffs_[i] = fitParams_[i + 5]; }
	intensityRatio_ = 4 / coeffs_[splineSizeZ_ / 2 * splineSizeXY_ + splineSizeY_ / 2 * splineSizeX_ + splineSizeX_ / 2];
}

//-------------------------------------------------------------------------------------------------
void SplineFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point)
{
	double offset    = 0,
	       intensity = 0,
	       x         = infos.x_,
	       y         = infos.y_,
	       z         = initZ_,
	       rss       = -1;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const size_t Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux <= width_ - roiHalf_ - 1 && roiHalf_ <= Uy && Uy <= height_ - roiHalf_ - 1) {
		// Création des variables
		double lambda = LAMBDA_0, minI, maxI, model, err, t1, t2;
		std::vector<double> hessian(nParams2_, 0.0),
		                    dudt(nParams_, 0.0), jacobian(nParams_, 0.0),
		                    data(n_);
		std::vector<size_t> coords(2 * n_);
		std::array<std::vector<double>, 3> p;

		// Initialization of data & t (A vérifier)
		preCompute(Ux, Uy, data, coords, minI, maxI);

		// Initialisation du/des P :
		// - Calcul du centre de masse de la région pour x et y
		CenterOfMass(data, roi_, x, y);

		initP(p, x, y, z, minI, maxI);
		p[1][1] *= intensityRatio_;         // Change l'intensité pour la position actuelle.
		p[2][1] = std::max(100.0, p[1][1]); // Change l'intenisté pour le saut maximal.

		// Calcul du patch 3D de 4x4x4 (64)
		const double Xc = p[1][2] - std::floor(p[1][2]),
		             Yc = p[1][3] - std::floor(p[1][3]),
		             Zc = p[1][4] - std::floor(p[1][4]);
		std::vector<double> f(64, 0.0), dxf(64, 0.0), dyf(64, 0.0), dzf(64, 0.0);
		Delta3D(Xc, Yc, Zc, f, dxf, dyf, dzf);

		// Calcul des dérivées et initialisation du model
		// offset utile, car la ROI et les coeffs ont peut-être des tailles différentes donc on recentre sur X et sur Y
		// Exemple : ROI de 7 splineX de 11 : La ROI 7x7 est centrée à partir du point 2 dans la spline
		const size_t Off = static_cast<size_t>(std::floor((static_cast<double>(splineSizeX_) + 1.0 - static_cast<double>(roi_)) / 2));
		for (size_t i = 0; i < n_; ++i) {
			derivativeSpline(coords[i * 2] + Off, coords[i * 2 + 1] + Off, static_cast<size_t>(std::floor(p[1][4])), f, dxf, dyf, dzf, p[1], dudt, model);
			updateModel(data[i], model, err, t1, t2, jacobian, hessian, dudt);
		}

		bool stop        = false;
		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			//...........................
			if (iteration >= N_ITER) { stop = true; }
		}
	}
	fill(infos, point, offset, intensity, x, y, z);
}

//-------------------------------------------------------------------------------------------------
void SplineFit::derivativeSpline(const size_t x, const size_t y, const size_t z, const std::vector<double>& f,
                                 const std::vector<double>& dxf, const std::vector<double>& dyf, const std::vector<double>& dzf,
                                 const std::vector<double>& p, std::vector<double>& dudt, double& model) const
{
	double temp = 0;
	std::ranges::fill(dudt, 0.0);
	const size_t Pos = z * splineSizeXY_ + y * splineSizeX_ + x;

	for (size_t i = 0; i < 64; i++) {
		const size_t Id = i * splineSizeXYZ_ + Pos;
		temp += f[i] * coeffs_[Id];
		dudt[2] += dxf[i] * coeffs_[Id];
		dudt[3] += dyf[i] * coeffs_[Id];
		dudt[4] += dzf[i] * coeffs_[Id];
	}

	dudt[0] = 1;          // offset
	dudt[1] = temp;       // Intensité
	dudt[2] *= -1 * p[1]; // x
	dudt[3] *= -1 * p[1]; // y
	dudt[4] *= p[1];      // z
	model = p[0] + p[1] * temp;
}

void SplineFit::updateModel(const double data, const double model, double& err, double& t1, double& t2,
                            std::vector<double>& jacobian, std::vector<double>& hessian, const std::vector<double>& dudt) const
{
	if (data > 0) { err = err + 2 * (model - data - data * log(model / data)); }
	else { err = err + 2 * model; }

	t1 = 1 - data / model;
	for (size_t i = 0; i < nParams_; ++i) { jacobian[i] += t1 * dudt[i]; }

	t2 = data / pow(model, 2);
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = i; j < nParams_; ++j) {
			hessian[i * nParams_ + j] += t2 * dudt[i] * dudt[j];
			hessian[j * nParams_ + i] = hessian[i * nParams_ + j];
		}
	}
}
}
