// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <numeric>
#include <ranges>
#include <thread>
#include <vector>

#include "ATrous.hpp"
#include "Segmentation2D.hpp"
#include "fit/GaussianFit.hpp"
#include "fit/SplineFit.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
Palm::Palm(const uint16_t* image, double* points, const size_t pointsListSize, const size_t height, const size_t width, const size_t planes,
           const double thresholdVal, const double watershedRatio, const size_t fit, const double* fitParams)
	: stack_(image), points_(points), height_(height), width_(width), planes_(planes), size_(height * width),
	  nPointsMaxPerPlane_((pointsListSize / N_LOCALIZATIONS) / planes),
	  thresholdVal_(thresholdVal), watershedRatio_(watershedRatio), fit_(fit), fitParams_(fitParams) { }

//-------------------------------------------------------------------------------------------------
int Palm::Processing() const
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique du nombre de threads.
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
double Palm::AutoThreshold()
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique du nombre de threads.
	std::vector<std::array<double, N_LOCALIZATIONS>> points;
	std::vector mask(size_, false); //  Creation du masque binaire (true = masque actif)
	const double HalfRoi = fitParams_[0] / 2;

	// Calcul de l'écart type initial
	double mu       = std::accumulate(stack_, stack_ + size_, 0.0) / static_cast<double>(size_);
	double variance = 0.0;
	for (size_t i = 0; i < size_; ++i) { variance += (stack_[i] - mu) * (stack_[i] - mu); }
	variance /= static_cast<double>(size_);
	thresholdVal_ = std::sqrt(variance);

	for (size_t i = 0; i < THRESHOLD_ITERATION; ++i) {
		processingPlane(0, points);
		// Mise à jour du masque
		std::fill(mask.begin(), mask.end(), false);
		for (auto& p : points) {
			const auto Xmin = static_cast<size_t>(std::max(0, static_cast<int>(p[4] - HalfRoi))),
			           Xmax = std::min(width_, static_cast<size_t>(p[4] + HalfRoi)),
			           Ymin = static_cast<size_t>(std::max(0, static_cast<int>(p[5] - HalfRoi))),
			           Ymax = std::min(height_, static_cast<size_t>(p[5] + HalfRoi));
			for (size_t y = Ymin; y < Ymax; ++y) { for (size_t x = Xmin; x < Xmax; ++x) { mask[y * width_ + x] = true; } }
		}
		// Extraction des pixels non masqués (hors ROI)
		std::vector<double> tmp;
		tmp.reserve(size_); // pré-allocation
		for (size_t j = 0; j < size_; ++j) {
			if (!mask[j]) {
				tmp.push_back(stack_[j]); // image_ : image 2D stockée en 1D
			}
		}

		// Si nombre de pixels restant insuffisant break
		if (tmp.empty()) { break; }

		// Calcul de l'écart type pour les pixels hors segmentation
		mu       = std::accumulate(tmp.begin(), tmp.end(), 0.0) / static_cast<double>(tmp.size());
		variance = 0.0;
		for (const double V : tmp) { variance += (V - mu) * (V - mu); }
		variance /= static_cast<double>(tmp.size());
		variance = std::sqrt(variance);
		if (std::fabs(thresholdVal_ - variance) <= std::numeric_limits<double>::epsilon()) { break; } // Aucun changement
		thresholdVal_ = variance;
	}

	return thresholdVal_;
}

//-------------------------------------------------------------------------------------------------
void Palm::processingPlane(const size_t plane, std::vector<std::array<double, N_LOCALIZATIONS>>& points) const
{
	const uint16_t* image = stack_ + plane * size_;

	//----- Décomposition en ondelette -----
	const ATrous Decomposition(image, height_, width_); // Objet pour la décomposition ondelettes (algorithme à trous).

	//----- Segmentation -----
	const Segmentation2D Segmentation(Decomposition, thresholdVal_, watershedRatio_); // Objet de segmentation intelligente 2D.

	//----- Fitting -----
	switch (fit_) {
		case 1: {
			GaussianFit fit(fit_mode::gaussian_2d_fixed, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(Segmentation, points);
			break;
		}

		case 2: { // Gaussian fit of 2D
			GaussianFit fit(fit_mode::gaussian_2d_sigma, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(Segmentation, points);
			break;
		}

		case 3: { // Gaussian fit of 3D Theta is constant
			GaussianFit fit(fit_mode::gaussian_3d_fixed, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(Segmentation, points);
			break;
		}

		case 4: { // Gaussian fit of 3D
			GaussianFit fit(fit_mode::gaussian_3d_theta, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(Segmentation, points);
			break;
		}

		case 5: { // Spline fit for Z on EMCCD sensor (maybe same call for sCMOS)
			SplineFit fit(fit_mode::spline_3d_emccd, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(Segmentation, points);
			break;
		}

		default: {
			Fitting fit(fit_mode::none, image, height_, width_, nPointsMaxPerPlane_, fitParams_);
			fit.Fit(Segmentation, points);
			break;
		}
	}

	// Tri par Y puis X (en cas peu probable de Y identique)
	std::ranges::sort(points, [](const auto& a, const auto& b) {
		if (std::fabs(a[5] - b[5]) > EPSILON) { return a[5] < b[5]; } // Tri sur Y (p[5])
		return a[4] < b[4];                                           // Puis sur X (p[4])
	});
}
} // namespace palm_tracer
