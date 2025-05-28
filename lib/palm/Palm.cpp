// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <numbers>
#include <thread>
#include <vector>

#include "ATrous.hpp"
#include "Segmentation2D.hpp"
#include "fit/GaussianFit.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
Palm::Palm(const unsigned short* image, double* points, const size_t pointsListSize, const size_t height, const size_t width, const size_t planes,
           const double thresholdVal, const double watershedRatio, const size_t fit, const double* fit_params)
	: stack_(image), points_(points), height_(height), width_(width), planes_(planes), size_(height * width),
	  nPointsMax_(pointsListSize / N_LOCALIZATIONS), nPointsMaxPerPlane_((pointsListSize / N_LOCALIZATIONS) / planes),
	  thresholdVal_(thresholdVal), watershedRatio_(watershedRatio), fit_(fit)
{
	if (1 <= fit_ && fit_ <= 4) { // Fit Gaussien
		roi_        = static_cast<size_t>(fit_params[0]);
		initSigmaX_ = fit_params[1];
		initSigmaY_ = 2 * fit_params[1];
		initTheta_  = fit_params[2];
	}
	//else if (5 <= fit_ && fit_ <= XX) {} //Fit Spline
	else { // Aucun Fit
		roi_        = 0;
		initSigmaX_ = 0;
		initSigmaY_ = 0;
		initTheta_  = 0;
	}
}

//-------------------------------------------------------------------------------------------------
int Palm::Processing() const
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique
	//omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max

	std::vector<std::vector<std::array<double, N_LOCALIZATIONS>>> allPoints(planes_);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(planes_); ++i) { processingPlane(i, allPoints[i]); }

	// Copy et remplissage Remplissage de l'index global, du plan et de l'index par plan
	size_t count = 0, global = 1;
	for (int i = 0; i < static_cast<int>(planes_); ++i) {
		size_t local = 1;
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
	Segmentation2D segmentation(decomposition, thresholdVal_, watershedRatio_);  // Objet de segmentation intelligente 2D.

	//----- Fitting -----
	//Fitting fit(fit_mode::none, image, height_, width_, nPointsMaxPerPlane_);
	//fit.Fit(segmentation, points);
	const auto& P = segmentation.Infos();
	points.resize(P.size());

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(std::min<size_t>(P.size(), nPointsMaxPerPlane_ - 1)); ++i) {
		double x         = P[i].x_;
		double y         = P[i].y_;
		double sigmaX    = initSigmaX_;
		double sigmaY    = initSigmaY_;
		double theta     = initTheta_;
		double intensity = 0, offset = 0, rss = -1;

		switch (fit_) {
			case 1: { // Gaussian fit of 2D Sigma is constant
				GaussianFit fit(fit_mode::gaussian_2d_fixed, image, height_, width_, nPointsMaxPerPlane_, initSigmaX_, initSigmaY_, initTheta_);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				sigmaX = initSigmaX_;
				sigmaY = initSigmaX_;
				theta  = 0;
				break;
			}

			case 2: { // Gaussian fit of 2D
				GaussianFit fit(fit_mode::gaussian_2d_sigma, image, height_, width_, nPointsMaxPerPlane_, initSigmaX_, initSigmaY_, initTheta_);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				sigmaY = sigmaX;
				theta  = 0;
				break;
			}

			case 3: { // Gaussian fit of 3D Theta is constant
				GaussianFit fit(fit_mode::gaussian_3d_fixed, image, height_, width_, nPointsMaxPerPlane_, initSigmaX_, initSigmaY_, initTheta_);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			case 4: { // Gaussian fit of 3D
				GaussianFit fit(fit_mode::gaussian_3d_theta, image, height_, width_, nPointsMaxPerPlane_, initSigmaX_, initSigmaY_, initTheta_);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			default: {
				sigmaX = sigmaY = theta = rss = 0;
				break;
			}
		}
		auto& p = points[i];
		p[0]    = 0;                              // Id
		p[1]    = static_cast<double>(plane + 1); // Plane
		p[2]    = 0;                              // Index
		p[3]    = -1;                             // Channel
		p[4]    = x;                              // X
		p[5]    = y;                              // Y
		p[6]    = 0;                              // Z

		// Integrated intensity
		if (fit_ == 0) { p[7] = P[i].intensity_; }
		else { p[7] = 2.0 * std::numbers::pi * intensity * sigmaX * sigmaY; }

		p[8]  = sigmaX;          // Sigma X
		p[9]  = sigmaY;          // Sigma Y
		p[10] = theta;           // Theta
		p[11] = rss;             // MSE Gaussian
		p[12] = -1;              // MSE Z
		p[13] = 0;               // Pair Distance
		p[14] = intensity;       // Intensity 0
		p[15] = offset;          // Intensity offset
		p[16] = P[i].intensity_; // Intensity
		p[17] = P[i].surf_;      // Surface

		// Circularity
		if (sigmaX > EPSILON && sigmaY > EPSILON) { p[18] = std::min(sigmaX, sigmaY) / std::max(sigmaX, sigmaY); }
		else { p[18] = 1; }
	}

	// Tri par Y puis X (en cas peu probable de Y identique)
	std::ranges::sort(points, [](const auto& a, const auto& b) {
		if (std::fabs(a[5] - b[5]) > EPSILON) { return a[5] < b[5]; } // Tri sur Y (p[4])
		return a[4] < b[4];                                           // Puis sur X (p[3])
	});
}
}
