///-------------------------------------------------------------------------------------------------
/// \file FrameProcessor.cpp
/// <summary>Implémentation du pipeline réutilisable de traitement d’une image.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#include "FrameProcessor.hpp"

#include <algorithm>
#include <compare>
#include <ranges>

#include "../fit/GaussianFit.hpp"
#include "../fit/SplineFit.hpp"

namespace palm_tracer::localization {
namespace col = localization_columns;

//-------------------------------------------------------------------------------------------------
FrameProcessor::FrameProcessor(const Vec2I& dim, const size_t maxPoints, const double watershedRatio,
                               const size_t fit, const double* fitParams)
	: dim_(dim), maxPoints_(maxPoints), watershedRatio_(watershedRatio), fit_(fit), fitParams_(fitParams),
	  decomposition_(dim), segmentation_(dim) { }

//-------------------------------------------------------------------------------------------------
void FrameProcessor::Process(const uint16_t* image, const double threshold, FrameLocalizations& points)
{
	decomposition_.Process(image);
	segmentation_.Process(decomposition_, threshold, watershedRatio_);
	fitRegions(image, points);

	// Tri par Y puis X (en cas peu probable de Y identique)
	std::ranges::sort(points, [](const auto& a, const auto& b) {
		// L'ordre faible reste valide même si un ajustement produit un NaN.
		if (const auto Order = std::weak_order(a[col::Y], b[col::Y]); Order != 0) { return Order < 0; }
		return std::weak_order(a[col::X], b[col::X]) < 0;
	});
}

//-------------------------------------------------------------------------------------------------
void FrameProcessor::fitRegions(const uint16_t* image, FrameLocalizations& points) const
{
	switch (fit_) {
		case 1: {
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_2d_fixed, image, dim_, maxPoints_, fitParams_);
			Gfit.FitAll(segmentation_, points);
			break;
		}

		case 2: { // Ajustement gaussien 2D
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_2d_sigma, image, dim_, maxPoints_, fitParams_);
			Gfit.FitAll(segmentation_, points);
			break;
		}

		case 3: { // Ajustement gaussien avec angle fixe
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_3d_fixed, image, dim_, maxPoints_, fitParams_);
			Gfit.FitAll(segmentation_, points);
			break;
		}

		case 4: { // Ajustement gaussien avec angle libre
			const fit::GaussianFit Gfit(fit::fit_mode::gaussian_3d_theta, image, dim_, maxPoints_, fitParams_);
			Gfit.FitAll(segmentation_, points);
			break;
		}

		case 5: { // Ajustement axial par spline pour un capteur EMCCD
			const fit::SplineFit Sfit(fit::fit_mode::spline_3d_emccd, image, dim_, maxPoints_, fitParams_);
			Sfit.FitAll(segmentation_, points);
			break;
		}

		default: {
			const fit::Fitting Sfit(fit::fit_mode::none, image, dim_, maxPoints_, fitParams_);
			Sfit.FitAll(segmentation_, points);
			break;
		}
	}
}
} // namespace palm_tracer::localization
