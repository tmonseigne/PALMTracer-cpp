// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.hpp"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <thread>

#include "ATrous.hpp"
#include "Segmentation2D.hpp"
#include "fit/GaussianFit.hpp"


//-------------------------------------------------------------------------------------------------
int Palm::Processing()
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique
	//omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max

	ATrous decomposition; // Objet pour la décomposition ondelettes (algorithme à trous).
	decomposition.ATrousDecomposition2D(image_, height_, width_);

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
	segmentation.ComputeStatistics(output, segmentMap, decomposition.wDecomp_[PLANE_ID], gaussFit_);
	//if (VOLUME_MIN > 0. || INTENSITY_MIN > 0.) { segmentation.FilterData(segmentMap, structFiltParam_.data(), structFiltParamUse_.data()); }
	segmentation.FilterData(segmentMap);

	const auto& P = segmentation.GetInfos();
	std::vector<std::array<double, N_SEGMENTS>> allPoints(P.size());

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(std::min<size_t>(P.size(), nPointsMaxPerPlane_ - 1)); ++i) {
		double x = P[i].moment_.y;
		double y = P[i].moment_.x;
		double sigmaX = initSigmaX_;
		double sigmaY = initSigmaY_;
		double theta = initTheta_;
		double intensity = 0, offset = 0, rss = -1;

		switch (gaussFit_) {
			case 1: { // Gaussian fit of 2D Sigma is constant
				GaussianFit fit(fit_mode::gauss_2d);
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				sigmaX = initSigmaX_;
				sigmaY = initSigmaX_;
				theta = 0;
				break;
			}

			case 2: { // Gaussian fit of 2D
				GaussianFit fit(fit_mode::gauss_2d_sigma);
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				sigmaY = sigmaX;
				theta = 0;
				break;
			}

			case 3: { // Gaussian fit of 3D Theta is constant
				GaussianFit fit(fit_mode::gauss_3d);
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			case 4: { // Gaussian fit of 3D
				GaussianFit fit(fit_mode::gauss_3d_theta);
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			default: {
				sigmaX = sigmaY = theta = rss = 0;
				break;
			}
		}
		auto& point = allPoints[i];
		point[0] = sigmaX;
		point[1] = sigmaY;
		point[2] = theta;
		point[3] = y;
		point[4] = x;
		point[5] = intensity;
		point[6] = offset;
		point[7] = rss;
		point[8] = P[i].intensity_;
		point[9] = P[i].surf_;
		point[10] = 0.0f; // Z coordinate
		point[11] = 0.0f; // Spectral or other
		point[12] = 0.0f; // Other
	}

	for (const auto& Point : allPoints) {
		std::ranges::copy(Point, points_ + currentPoint_);
		currentPoint_ += N_SEGMENTS;
	}
	// Blank Line
	std::fill_n(points_ + currentPoint_, N_SEGMENTS, -1);
	currentPoint_ += N_SEGMENTS;

	return static_cast<int>(currentPoint_);
}
