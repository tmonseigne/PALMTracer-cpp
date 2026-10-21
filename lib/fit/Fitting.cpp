#include "Fitting.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <iostream>
#include <numbers>
#include <vector>

namespace palm_tracer::fit {
namespace col = localization_columns;

//-------------------------------------------------------------------------------------------------
void Fitting::initialize()
{
	// Détermination du nombre de paramètres selon le modèle choisi
	switch (mode_) {
		case fit_mode::none:
		case fit_mode::gaussian_2d_fixed: nParams_ = 4; // offset, amplitude, x, y
			break;
		case fit_mode::gaussian_2d_sigma: nParams_ = 5; // offset, amplitude, x, y, sigma
			break;
		case fit_mode::gaussian_3d_fixed: nParams_ = 6; // offset, amplitude, x, y, sigmaX, sigmaY
			break;
		case fit_mode::gaussian_3d_theta: nParams_ = 7; // offset, amplitude, x, y, sigmaX, sigmaY, theta
			break;
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: nParams_ = 5; // offset, amplitude, x, y, z
			break;
	}
	nParams2_ = nParams_ * nParams_;

	// Paramètres utilisateur → membres
	parseParams();

	// ROI : forcer une taille impaire
	if (roi_ % 2 == 0) { roi_--; }
	roiHalf_ = roi_ / 2;
	n_       = roi_ * roi_;
}

// ==================== Fit Computes ====================
//-------------------------------------------------------------------------------------------------
void Fitting::FitAll(const localization::Segmentation2D& segmentation, std::vector<std::array<double, col::Count>>& points) const
{
	const auto& P      = segmentation.Infos();
	const size_t Count = std::min(P.size(), maxPoints_);
	if (P.size() > maxPoints_) {
		#pragma omp critical(FittingCapacityWarning)
		{
			std::cerr << "[WARNING] Density is too high; some points are missed (output capacity exceeded).\n";
		}
	}
	points.resize(Count);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(Count); ++i) { Fit(P[i], points[i]); }
}

//-------------------------------------------------------------------------------------------------
void Fitting::fill(const localization::SegmentationInfo& infos, std::array<double, col::Count>& point, const double x, const double y, const double z,
                   const double offset, const double intensity, const double sigmaX, const double sigmaY, const double theta,
                   const double rssXY, const double rssZ) const
{
	// Champs d'identité (non gérés ici)
	point[col::Id]      = 0;  // Id
	point[col::Plane]   = 0;  // Plane
	point[col::Index]   = 0;  // Index
	point[col::Channel] = -1; // Channel

	// Position (pixels ou unités adaptées pour z)
	point[col::X] = x; // X
	point[col::Y] = y; // Y
	point[col::Z] = z; // Z

	// Intensité intégrée selon le modèle
	switch (mode_) {
		case fit_mode::none: point[col::Intensity] = infos.intensity_;
			break;
		case fit_mode::gaussian_2d_fixed:
		case fit_mode::gaussian_2d_sigma:
		case fit_mode::gaussian_3d_fixed:
		case fit_mode::gaussian_3d_theta:
			// I_int ≈ 2π × amplitude × sigmaX × sigmaY (pixels)
			if (intensity > EPSILON && sigmaX > EPSILON && sigmaY > EPSILON) { point[col::Intensity] = 2.0 * std::numbers::pi * intensity * sigmaX * sigmaY; }
			else { point[col::Intensity] = -1; }
			break;
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: point[col::Intensity] = intensity;
			break;
	}

	point[col::SigmaX]       = sigmaX;           // Sigma X
	point[col::SigmaY]       = sigmaY;           // Sigma Y
	point[col::Theta]        = theta;            // Theta en radians pour le modèle gaussien.
	point[col::QualityXY]    = rssXY;            // Indicateur XY ; colonne historique MSE XY.
	point[col::QualityZ]     = rssZ;             // Indicateur axial ; colonne historique MSE Z.
	point[col::Amplitude]    = intensity;        // Amplitude
	point[col::Offset]       = offset;           // Offset
	point[col::SegIntensity] = infos.intensity_; // Intensité de la segmentation
	point[col::Surface]      = infos.surf_;      // Surface de la segmentation

	// Circularité = min(sigmaX, sigmaY) / max(sigmaX, sigmaY) si définis, 1 sinon
	if (sigmaX > EPSILON && sigmaY > EPSILON) { point[col::Circularity] = std::min(sigmaX, sigmaY) / std::max(sigmaX, sigmaY); }
	else { point[col::Circularity] = 1; }
}

//-------------------------------------------------------------------------------------------------
void Fitting::preCompute(const size_t x, const size_t y, std::vector<Pixel>& grid, double& min, double& max) const
{
	min = std::numeric_limits<double>::max();    // Plus grand double possible (DBL_MAX)
	max = std::numeric_limits<double>::lowest(); // Plus petit double possible (-DBL_MAX)

	const int Half = static_cast<int>(roiHalf_);
	for (int i = -Half; i <= Half; ++i) {
		const auto Idx = Half + ((i + Half) * roi_);
		for (int j = -Half; j <= Half; ++j) {
			grid[j + Idx] = {.x = x + j, .y = y + i, .v = static_cast<double>(image_[x + j + ((y + i) * dim_.x)])};
			min           = std::min(min, grid[j + Idx].v);
			max           = std::max(max, grid[j + Idx].v);
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Fitting::initP(std::array<std::vector<double>, 3>& p, const double x, const double y, const double z, const double minI, const double maxI) const
{
	constexpr double MinI = 0.0, MaxI = 65535.0,
	                 MinS = 0.1, MaxS = 20.0;

	switch (mode_) {
		case fit_mode::none:
		case fit_mode::gaussian_2d_fixed: {
			// [offset, amplitude, x, y]
			p[0] = {MinI, MinI, x - 1, y - 1};
			p[1] = {minI, maxI - minI, x, y};
			p[2] = {MaxI, MaxI, x + 1, y + 1};
			break;
		}
		case fit_mode::gaussian_2d_sigma: {
			// [offset, amplitude, x, y, sigma]
			p[0] = {MinI, MinI, x - 1, y - 1, MinS};
			p[1] = {minI, maxI - minI, x, y, fitParams_[1]};
			p[2] = {MaxI, MaxI, x + 1, y + 1, MaxS};
			break;
		}
		case fit_mode::gaussian_3d_fixed: {
			// [offset, amplitude, x, y, sigmaX, sigmaY]
			p[0] = {MinI, MinI, x - 1, y - 1, MinS, MinS};
			p[1] = {minI, maxI - minI, x, y, std::numbers::sqrt2_v<double> * fitParams_[1], std::numbers::sqrt2_v<double> * fitParams_[2]};
			p[2] = {MaxI, MaxI, x + 1, y + 1, MaxS, MaxS};
			break;
		}
		case fit_mode::gaussian_3d_theta: {
			// [offset, amplitude, x, y, sigmaX, sigmaY, theta]
			p[0] = {MinI, MinI, x - 1, y - 1, MinS, MinS, -150};
			p[1] = {minI, maxI - minI, x, y, std::numbers::sqrt2_v<double> * fitParams_[1], std::numbers::sqrt2_v<double> * fitParams_[2], fitParams_[3]};
			p[2] = {MaxI, MaxI, x + 1, y + 1, MaxS, MaxS, 150};
			break;
		}
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: {
			// [offset, amplitude, x, y, z]
			const double RoiC = (fitParams_[0] - 1) / 2.0,
			             RoiD = fitParams_[0] / 4.0;
			const double Dmin = RoiC - RoiD,
			             Dmax = RoiC + RoiD;

			p[0] = {-0.01, -1, Dmin, Dmin, 0};
			p[1] = {minI, maxI - minI, x, y, z};
			p[2] = {0.01, 1, Dmax, Dmax, 1}; // Le max pour Z est borné par la taille des PSF spline (voir paramètres).
			break;
		}
	}
}
} // namespace palm_tracer::fit
