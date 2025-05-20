// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.hpp"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <ranges>
#include <thread>

#include "fit/GaussianFit.hpp"


//-------------------------------------------------------------------------------------------------
Palm::Palm(const unsigned short* image, double* points, const size_t pointsListSize, const size_t height, const size_t width,
           const double thresholdVal, const double watershedRatio, const size_t gaussFit, const double sigma, const double theta, const size_t roi)
	: image_(image), points_(points), height_(height), width_(width), size_(static_cast<size_t>(height) * width),
	  nPointsMaxPerPlane_(pointsListSize / N_SEGMENTS), currentPoint_(0), roi_(roi), thresholdVal_(thresholdVal),
	  watershedRatio_(watershedRatio), gaussFit_(gaussFit), initSigmaX_(sigma), initSigmaY_(2 * sigma), initTheta_(theta),
	  segmentPub_(height, width)
{
	operators_.fill(0);
	segPub_.resize(size_);
	outputPub_.resize(size_);
}


//-------------------------------------------------------------------------------------------------
int Palm::Processing()
{
	omp_set_dynamic(1); // Autorise l’adaptation automatique
	//omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max

	aTrousPub_.ATrousDecomposition2D(image_, height_, width_);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(size_); ++i) {
		if (aTrousPub_.wDecomp_[PLANE_ID][i] > thresholdVal_) {
			outputPub_[i] = aTrousPub_.wDecomp_[PLANE_ID][i];
			segPub_[i] = 1;
		}
		else {
			segPub_[i] = 0;
			outputPub_[i] = 0;
		}
	}

	segmentPub_.WatershedIntelligent(outputPub_, segPub_, watershedRatio_);
	segmentPub_.ComputeStatistics(outputPub_, segPub_, aTrousPub_.wDecomp_[PLANE_ID], gaussFit_);
	//if (VOLUME_MIN > 0. || INTENSITY_MIN > 0.) { segmentPub_.FilterData(segPub_, structFiltParam_.data(), structFiltParamUse_.data()); }
	segmentPub_.FilterData(segPub_);

	const auto& P = segmentPub_.GetInfos();
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
