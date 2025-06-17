// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <thread>
#include <vector>

#include "ATrous.hpp"
#include "Segmentation2D.hpp"
#include "fit/GaussianFit.hpp"
#include "fit/SplineFit.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
Palm::Palm(const unsigned short* image, double* points, const size_t pointsListSize, const size_t height, const size_t width, const size_t planes,
           const double thresholdVal, const double watershedRatio, const size_t fit, const double* fitParams)
	: stack_(image), points_(points), height_(height), width_(width), planes_(planes), size_(height * width),
	  nPointsMaxPerPlane_((pointsListSize / N_LOCALIZATIONS) / planes),
	  thresholdVal_(thresholdVal), watershedRatio_(watershedRatio), fit_(fit), fitParams_(fitParams) { }

//-------------------------------------------------------------------------------------------------
int Palm::Processing() const
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique du nombre de thread
	//omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max

	std::vector<std::vector<std::array<double, N_LOCALIZATIONS>>> allPoints(planes_);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(planes_); ++i) { processingPlane(i, allPoints[i]); }

	// Copy et remplissage Remplissage de l'index global, du plan et de l'index par plan
	size_t count = 0, global = 1;
	for (int i = 0; i < static_cast<int>(planes_); ++i) {
		size_t local       = 1;
		const double Plane = i + 1;
		for (const auto& Point : allPoints[i]) {
			std::ranges::copy(Point, points_ + count);
			points_[count]     = static_cast<double>(global);
			points_[count + 1] = Plane;
			points_[count + 2] = static_cast<double>(local);
			++global;
			++local;
			count += N_LOCALIZATIONS;
		}
	}

	return static_cast<int>(count);
}

//-------------------------------------------------------------------------------------------------
void Palm::processingPlane(const size_t plane, std::vector<std::array<double, N_LOCALIZATIONS>>& points) const
{
	const unsigned short* image = stack_ + (plane * size_);

	//----- Décomposition en ondelette -----
	ATrous decomposition(image, height_, width_); // Objet pour la décomposition ondelettes (algorithme à trous).

	//----- Segmentation -----
	Segmentation2D segmentation(decomposition, thresholdVal_, watershedRatio_); // Objet de segmentation intelligente 2D.

	//----- Fitting -----
	switch (fit_) {
		case 1: {
			GaussianFit fit(fit_mode::gaussian_2d_fixed, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(segmentation, points);
			break;
		}

		case 2: { // Gaussian fit of 2D
			GaussianFit fit(fit_mode::gaussian_2d_sigma, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(segmentation, points);
			break;
		}

		case 3: { // Gaussian fit of 3D Theta is constant
			GaussianFit fit(fit_mode::gaussian_3d_fixed, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(segmentation, points);
			break;
		}

		case 4: { // Gaussian fit of 3D
			GaussianFit fit(fit_mode::gaussian_3d_theta, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(segmentation, points);
			break;
		}

		case 5: { // Spline fit for Z on EMCCD sensor (maybe same call for sCMOS)
			SplineFit fit(fit_mode::spline_3d_emccd, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(segmentation, points);
			break;
		}

		default: {
			Fitting fit(fit_mode::none, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(segmentation, points);
			break;
		}
	}

	// Tri par Y puis X (en cas peu probable de Y identique)
	std::ranges::sort(points, [](const auto& a, const auto& b) {
		if (std::fabs(a[5] - b[5]) > EPSILON) { return a[5] < b[5]; } // Tri sur Y (p[5])
		return a[4] < b[4];                                           // Puis sur X (p[4])
	});
}
}
