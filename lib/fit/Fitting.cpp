#include "Fitting.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <numbers>
#include <vector>

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Fitting::initialize()
{
	switch (mode_) {
		case fit_mode::none:
		case fit_mode::gaussian_2d_fixed: nParams_ = 4; // 4 paramètres : offset, amplitude, x, y
			break;
		case fit_mode::gaussian_2d_sigma: nParams_ = 5; // 5 : offset, amplitude, x, y, sigma
			break;
		case fit_mode::gaussian_3d_fixed: nParams_ = 6; // 6 : offset, amplitude, x, y, sigmaX, sigmaY
			break;
		case fit_mode::gaussian_3d_theta: nParams_ = 7; // 7 : offset, amplitude, x, y, sigmaX, sigmaY, theta
			break;
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: nParams_ = 5; // 5 : offset, amplitude, x, y, z
			break;
	}
	nParams2_ = nParams_ * nParams_;
	parseParams();
	if (roi_ % 2 == 0) { roi_--; }
	roiHalf_ = roi_ / 2;
	n_       = roi_ * roi_;
}

// ==================== Fit Computes ====================
//-------------------------------------------------------------------------------------------------
void Fitting::Fit(const Segmentation2D& segmentation, std::vector<std::array<double, N_LOCALIZATIONS>>& points) const
{
	const auto& P = segmentation.Infos();
	points.resize(P.size());

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(std::min<size_t>(P.size(), maxPoints_ - 1)); ++i) { fit(P[i], points[i]); }
}

//-------------------------------------------------------------------------------------------------
void Fitting::fill(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point, const double offset, const double intensity,
                   const double x, const double y, const double z, const double sigmaX, const double sigmaY, const double theta, const double rss,
                   const double rssZ) const
{
	point[0] = 0;  // Id
	point[1] = 0;  // Plane
	point[2] = 0;  // Index
	point[3] = -1; // Channel
	point[4] = x;  // X
	point[5] = y;  // Y
	point[6] = z;  // Z
	// Integrated intensity
	switch (mode_) {
		case fit_mode::none: point[7] = infos.intensity_;
			break;
		case fit_mode::gaussian_2d_fixed:
		case fit_mode::gaussian_2d_sigma:
		case fit_mode::gaussian_3d_fixed:
		case fit_mode::gaussian_3d_theta: point[7] = 2.0 * std::numbers::pi * intensity * sigmaX * sigmaY;
			break;
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: point[7] = intensity;

			break;
	}

	point[8]  = sigmaX;           // Sigma X
	point[9]  = sigmaY;           // Sigma Y
	point[10] = theta;            // Theta
	point[11] = rss;              // MSE X, Y
	point[12] = rssZ;             // MSE Z
	point[13] = intensity;        // Intensity 0
	point[14] = offset;           // Intensity offset
	point[15] = infos.intensity_; // Intensity
	point[16] = infos.surf_;      // Surface

	// Circularity
	if (sigmaX > EPSILON && sigmaY > EPSILON) { point[17] = std::min(sigmaX, sigmaY) / std::max(sigmaX, sigmaY); }
	else { point[17] = 1; }
}

//-------------------------------------------------------------------------------------------------
void Fitting::preCompute(const size_t x, const size_t y, std::vector<Pixel>& grid, double& min, double& max) const
{
	min = std::numeric_limits<double>::max();    // Plus grand double possible (DBL_MAX)
	max = std::numeric_limits<double>::lowest(); // Plus petit double possible (-DBL_MAX)

	const int Half = static_cast<int>(roiHalf_);
	for (int i = -Half; i <= Half; ++i) {
		const auto Idx = Half + (i + Half) * roi_;
		for (int j = -Half; j <= Half; ++j) {
			grid[j + Idx] = {.x = x + j, .y = y + i, .v = static_cast<double>(image_[x + j + (y + i) * dim_.x])};
			min           = std::min(min, grid[j + Idx].v);
			max           = std::max(max, grid[j + Idx].v);
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Fitting::initP(std::array<std::vector<double>, 3>& p, const double x, const double y, const double z, const double minI, const double maxI) const
{
	constexpr double MinI = 0.0, MaxI = 65535,
	                 MinS = 0.1, MaxS = 20;
	switch (mode_) {
		case fit_mode::none:
		case fit_mode::gaussian_2d_fixed: {
			p[0] = {MinI, MinI, x - 1, y - 1}; // minimum théorique
			p[1] = {minI, maxI - minI, x, y};  // Valeur initiale
			p[2] = {MaxI, MaxI, x + 1, y + 1}; // maximum théorique
			break;
		}
		case fit_mode::gaussian_2d_sigma: {
			p[0] = {MinI, MinI, x - 1, y - 1, MinS};
			p[1] = {minI, maxI - minI, x, y, fitParams_[1]};
			p[2] = {MaxI, MaxI, x + 1, y + 1, MaxS};
			break;
		}
		case fit_mode::gaussian_3d_fixed: {
			p[0] = {MinI, MinI, x - 1, y - 1, MinS, MinS};
			p[1] = {minI, maxI - minI, x, y, std::numbers::sqrt2_v<double> * fitParams_[1], std::numbers::sqrt2_v<double> * fitParams_[2]};
			p[2] = {MaxI, MaxI, x + 1, y + 1, MaxS, MaxS};
			break;
		}
		case fit_mode::gaussian_3d_theta: {
			p[0] = {MinI, MinI, x - 1, y - 1, MinS, MinS, -150};
			p[1] = {minI, maxI - minI, x, y, std::numbers::sqrt2_v<double> * fitParams_[1], std::numbers::sqrt2_v<double> * fitParams_[2], fitParams_[3]};
			p[2] = {MaxI, MaxI, x + 1, y + 1, MaxS, MaxS, 150};
			break;
		}
		case fit_mode::spline_3d_emccd:
		case fit_mode::spline_3d_scmos: {
			const double RoiC = (fitParams_[0] - 1) / 2.0,
			             RoiD = fitParams_[0] / 4.0;
			const double Dmin = RoiC - RoiD,
			             Dmax = RoiC + RoiD;

			p[0] = {MinI, MinI, Dmin, Dmin, 0};
			p[1] = {minI, maxI - minI, x, y, z};
			p[2] = {MaxI, MaxI, Dmax, Dmax, 1}; // le maximum pour Z est défini dans les paramètres par splineSizeZ_.
			break;
		}
	}
}
} // namespace palm_tracer
