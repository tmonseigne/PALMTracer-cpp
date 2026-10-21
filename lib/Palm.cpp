#include "Palm.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <execution>
#include <numeric>
#include <ranges>
#include <thread>
#include <vector>

#include "fit/GaussianFit.hpp"
#include "fit/SplineFit.hpp"
#include "localization/ATrous.hpp"
#include "localization/Segmentation2D.hpp"


namespace palm_tracer {
namespace col = localization_columns;

//-------------------------------------------------------------------------------------------------
Palm::Palm(const uint16_t* image, double* points, const size_t pointsListSize, const Vec3I& dim,
           const double thresholdVal, const double watershedRatio, const size_t fit, const double* fitParams)
	: stack_(image), points_(points), dim_(dim), size_(dim.GetPlaneSize()),
	  nPointsMaxPerPlane_(dim.z == 0 ? 0 : (pointsListSize / col::Count) / dim.z),
	  thresholdVal_(thresholdVal), watershedRatio_(watershedRatio), fit_(fit), fitParams_(fitParams) { }

//-------------------------------------------------------------------------------------------------
int Palm::Processing() const
{
	if (nPointsMaxPerPlane_ == 0) { return 0; }
	std::vector<std::vector<std::array<double, col::Count>>> allPoints(dim_.z);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(dim_.z); ++i) { processingPlane(i, allPoints[i]); }

	// Copie par plan et attribution de l'identifiant global, du plan et de l'index local.
	size_t count = 0, global = 1;
	for (int i = 0; i < static_cast<int>(dim_.z); ++i) {
		size_t local       = 1;
		const double Plane = i + 1;
		for (const auto& Point : allPoints[i]) {
			std::ranges::copy(Point, points_ + count);
			points_[count + col::Id]    = static_cast<double>(global);
			points_[count + col::Plane] = Plane;
			points_[count + col::Index] = static_cast<double>(local);
			++global;
			++local;
			count += col::Count;
		}
	}

	return static_cast<int>(count);
}

//-------------------------------------------------------------------------------------------------
double Palm::AutoThreshold()
{
	std::vector<std::array<double, col::Count>> points;
	std::vector mask(size_, false); //  Création du masque binaire (true = masque actif)
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
		for (const auto& P : points) {
			const auto Xmin = static_cast<size_t>(std::max(0, static_cast<int>(P[col::X] - HalfRoi))),
			           Xmax = std::min(dim_.x, static_cast<size_t>(P[col::X] + HalfRoi)),
			           Ymin = static_cast<size_t>(std::max(0, static_cast<int>(P[col::Y] - HalfRoi))),
			           Ymax = std::min(dim_.y, static_cast<size_t>(P[col::Y] + HalfRoi));
			for (size_t y = Ymin; y < Ymax; ++y) { for (size_t x = Xmin; x < Xmax; ++x) { mask[(y * dim_.x) + x] = true; } }
		}
		// Deux passages dans le même ordre, sans recopier les pixels hors masque.
		size_t count = 0;
		double sum   = 0;
		for (size_t j = 0; j < size_; ++j) {
			if (!mask[j]) {
				sum += stack_[j];
				++count;
			}
		}

		// Si nombre de pixels restant insuffisant break
		if (count == 0) { break; }

		// Calcul de l'écart type pour les pixels hors segmentation
		mu       = sum / static_cast<double>(count);
		variance = 0.0;
		for (size_t j = 0; j < size_; ++j) { if (!mask[j]) { variance += (stack_[j] - mu) * (stack_[j] - mu); } }
		variance /= static_cast<double>(count);
		variance = std::sqrt(variance);
		if (std::abs(thresholdVal_ - variance) <= std::numeric_limits<double>::epsilon()) { break; } // Aucun changement
		thresholdVal_ = variance;
	}

	return thresholdVal_;
}

//-------------------------------------------------------------------------------------------------
void Palm::processingPlane(const size_t plane, std::vector<std::array<double, col::Count>>& points) const
{
	const uint16_t* image = stack_ + (plane * size_);
	const Vec2I Dim       = {.x = dim_.x, .y = dim_.y};

	//----- Décomposition en ondelette -----
	const localization::ATrous Decomposition(image, Dim); // Objet pour la décomposition ondelettes (algorithme à trous).

	//----- Segmentation -----
	const localization::Segmentation2D Segmentation(Decomposition, thresholdVal_, watershedRatio_); // Objet de segmentation intelligente 2D.

	//----- Fitting -----
	switch (fit_) {
		case 1: {
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_2d_fixed, image, Dim, nPointsMaxPerPlane_, fitParams_);
			Gfit.FitAll(Segmentation, points);
			break;
		}

		case 2: { // Gaussian fit of 2D
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_2d_sigma, image, Dim, nPointsMaxPerPlane_, fitParams_);
			Gfit.FitAll(Segmentation, points);
			break;
		}

		case 3: { // Gaussian fit of 3D Theta is constant
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_3d_fixed, image, Dim, nPointsMaxPerPlane_, fitParams_);
			Gfit.FitAll(Segmentation, points);
			break;
		}

		case 4: { // Gaussian fit of 3D
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_3d_theta, image, Dim, nPointsMaxPerPlane_, fitParams_);
			Gfit.FitAll(Segmentation, points);
			break;
		}

		case 5: { // Spline fit for Z on EMCCD sensor (maybe same call for sCMOS)
			const fit::SplineFit Sfit(fit::fit_mode::spline_3d_emccd, image, Dim, nPointsMaxPerPlane_, fitParams_);
			Sfit.FitAll(Segmentation, points);
			break;
		}

		default: {
			const fit::Fitting Sfit(fit::fit_mode::none, image, Dim, nPointsMaxPerPlane_, fitParams_);
			Sfit.FitAll(Segmentation, points);
			break;
		}
	}

	// Tri par Y puis X (en cas peu probable de Y identique)
	std::ranges::sort(points, [](const auto& a, const auto& b) {
		// L'ordre faible reste valide même si un ajustement produit un NaN.
		if (const auto Order = std::weak_order(a[col::Y], b[col::Y]); Order != 0) { return Order < 0; }
		return std::weak_order(a[col::X], b[col::X]) < 0;
	});
}
} // namespace palm_tracer
