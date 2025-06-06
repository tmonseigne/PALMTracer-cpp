///-------------------------------------------------------------------------------------------------
///
/// \file SplineFit.hpp
/// \brief Classe d'ajustement par spline.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <vector>

#include "Fitting.hpp"

namespace palm_tracer {
/// <summary>Classe d'ajustement gaussien.</summary>
class SplineFit final : public Fitting {
public:
	explicit SplineFit(const fit_mode mode, const unsigned short* image, const size_t height, const size_t width, const size_t maxPoints,
	                   const size_t roi, const double sigmaX, const double sigmaY, const double theta)
		: Fitting(mode, image, height, width, maxPoints, roi, sigmaX, sigmaY, theta) { }

protected:
	/// <summary>Exécute un ajustement de modèle spline sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) override;
};
}
