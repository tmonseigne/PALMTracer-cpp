#include "Fitting.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <numbers>
#include <vector>

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
Fitting::Fitting(const fit_mode mode, const unsigned short* image, const size_t height, const size_t width, const size_t maxPoints,
                 const size_t roi, const double sigmaX, const double sigmaY, const double theta)
	: mode_(mode), image_(image), n_(0), nParams_(4), height_(height), width_(width), maxPoints_(maxPoints),
	  roi_(roi), sigmaX_(sigmaX), sigmaY_(sigmaY), theta_(theta)
{
	switch (mode_) {
		case fit_mode::none:
		case fit_mode::gaussian_2d_fixed: nParams_ = 4; // 4 paramètres : x, y, amplitude, offset
			break;
		case fit_mode::gaussian_2d_sigma: nParams_ = 5; // 5 : x, y, amplitude, sigma, offset
			break;
		case fit_mode::gaussian_3d_fixed: nParams_ = 6; // 6 : x, y, amplitude, sigmaX, sigmaY, offset
			break;
		case fit_mode::gaussian_3d_theta: nParams_ = 7; // 7 : x, y, amplitude, sigmaX, sigmaY, theta, offset
			break;
		case fit_mode::spline_2d: break;
	}
	nParams2_ = nParams_ * nParams_;
	if (roi_ % 2 == 0) { roi_--; }
	roiHalf_ = roi_ / 2;
	n_       = roi_ * roi_;
}

// ==================== Fit Computes ====================
//-------------------------------------------------------------------------------------------------
void Fitting::Fit(const Segmentation2D& segmentation, std::vector<std::array<double, N_LOCALIZATIONS>>& points)
{
	const auto& P = segmentation.Infos();
	points.resize(P.size());

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(std::min<size_t>(P.size(), maxPoints_ - 1)); ++i) { fit(P[i], points[i]); }
}

//-------------------------------------------------------------------------------------------------
void Fitting::fill(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point, const double x, const double y, const double intensity,
                   const double sigmaX, const double sigmaY, const double theta, const double rss, const double offset) const
{
	point[0] = 0;  // Id
	point[1] = 0;  // Plane
	point[2] = 0;  // Index
	point[3] = -1; // Channel
	point[4] = x;  // X
	point[5] = y;  // Y
	point[6] = 0;  // Z

	// Integrated intensity
	if (mode_ == fit_mode::none) { point[7] = infos.intensity_; }
	else { point[7] = 2.0 * std::numbers::pi * intensity * sigmaX * sigmaY; }

	point[8]  = sigmaX;           // Sigma X
	point[9]  = sigmaY;           // Sigma Y
	point[10] = theta;            // Theta
	point[11] = rss;              // MSE X, Y
	point[12] = -1;               // MSE Z
	point[13] = 0;                // Pair Distance
	point[14] = intensity;        // Intensity 0
	point[15] = offset;           // Intensity offset
	point[16] = infos.intensity_; // Intensity
	point[17] = infos.surf_;      // Surface

	// Circularity
	if (sigmaX > EPSILON && sigmaY > EPSILON) { point[18] = std::min(sigmaX, sigmaY) / std::max(sigmaX, sigmaY); }
	else { point[18] = 1; }
}
}
