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


//-------------------------------------------------------------------------------------------------
int Palm::Processing() const
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique
	//omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max

	std::vector<std::vector<std::array<double, N_SEGMENTS>>> allPoints(planes_);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(planes_); ++i) { processingPlane(i, allPoints[i]); }

	// Copy et remplissage Remplissage de l'index global et de l'index par plan
	size_t count = 0, global = 1;
	for (int i = 0; i < static_cast<int>(planes_); ++i) {
		size_t local = 1;
		for (const auto& Point : allPoints[i]) {
			std::ranges::copy(Point, points_ + count);
			points_[count] = static_cast<double>(global);
			points_[count + 2] = static_cast<double>(local);
			++global;
			++local;
			count += N_SEGMENTS;
		}
	}

	return static_cast<int>(count);
}

//-------------------------------------------------------------------------------------------------
void Palm::processingPlane(const size_t plane, std::vector<std::array<double, N_SEGMENTS>>& points) const
{
	const unsigned short* image = image_ + (plane * size_);

	//----- Décomposition en ondelette -----
	ATrous decomposition; // Objet pour la décomposition ondelettes (algorithme à trous).
	decomposition.ATrousDecomposition2D(image, height_, width_);

	std::vector<unsigned short> segmentMap(size_); // Image segmentée (carte binaire des objets détectés).
	std::vector<double> output(size_);             // Carte des valeurs utilisées pour la segmentation (après seuillage).
	Segmentation2D segmentation(height_, width_);  // Objet de segmentation intelligente 2D.

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(size_); ++i) {
		if (decomposition.wDecomp_[PLANE_ID][i] > thresholdVal_) {
			output[i] = decomposition.wDecomp_[PLANE_ID][i];
			segmentMap[i] = 1;
		}
		else {
			segmentMap[i] = 0;
			output[i] = 0;
		}
	}

	segmentation.WatershedIntelligent(output, segmentMap, watershedRatio_);
	segmentation.ComputeStatistics(output, segmentMap, decomposition.wDecomp_[PLANE_ID], fit_);
	segmentation.FilterData(segmentMap);

	const auto& P = segmentation.GetInfos();
	points.resize(P.size());

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(std::min<size_t>(P.size(), nPointsMaxPerPlane_ - 1)); ++i) {
		double x = P[i].moment_.y;
		double y = P[i].moment_.x;
		double sigmaX = initSigmaX_;
		double sigmaY = initSigmaY_;
		double theta = initTheta_;
		double intensity = 0, offset = 0, rss = -1;

		switch (fit_) {
			case 1: { // Gaussian fit of 2D Sigma is constant
				GaussianFit fit(fit_mode::gauss_2d);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				sigmaX = initSigmaX_;
				sigmaY = initSigmaX_;
				theta = 0;
				break;
			}

			case 2: { // Gaussian fit of 2D
				GaussianFit fit(fit_mode::gauss_2d_sigma);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				sigmaY = sigmaX;
				theta = 0;
				break;
			}

			case 3: { // Gaussian fit of 3D Theta is constant
				GaussianFit fit(fit_mode::gauss_3d);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			case 4: { // Gaussian fit of 3D
				GaussianFit fit(fit_mode::gauss_3d_theta);
				fit.Run(image, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			default: {
				sigmaX = sigmaY = theta = rss = 0;
				break;
			}
		}
		auto& p = points[i];
		p[0] = 0;                              // Id
		p[1] = static_cast<double>(plane + 1); // Plane
		p[2] = 0;                              // Index
		p[3] = x;                              // X
		p[4] = y;                              // Y
		p[5] = 0;                              // Z

		// Integrated intensity
		if (fit_ == 0) { p[6] = P[i].intensity_; }
		else { p[6] = 2.0 * std::numbers::pi * intensity * sigmaX * sigmaY; }

		p[7] = sigmaX;           // Sigma X
		p[8] = sigmaY;           // Sigma Y
		p[9] = theta;            // Theta
		p[10] = intensity;       // Intensity 0
		p[11] = offset;          // Intensity offset
		p[12] = P[i].intensity_; // Intensity
		p[13] = P[i].surf_;      // Surface
		p[14] = rss;             // MSE Gaussian

		// Circularity
		if (sigmaX > EPSILON && sigmaY > EPSILON) { p[15] = std::min(sigmaX, sigmaY) / std::max(sigmaX, sigmaY); }
		else { p[15] = 1; }
	}

	// Tri par Y puis X (en cas peu probable de Y identique)
	std::ranges::sort(points, [](const auto& a, const auto& b) {
		if (std::fabs(a[4] - b[4]) > EPSILON) { return a[4] < b[4]; } // Tri sur Y (p[4])
		return a[3] < b[3];                                           // Puis sur X (p[3])
	});
}
