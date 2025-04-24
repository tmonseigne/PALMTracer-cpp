// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.h"
#include <omp.h>
#include <thread>

#include "fit/GaussianFit2D.h"
#include "fit/GaussianFit2D_Sigma.h"
#include "fit/GaussianFit3D.h"
#include "fit/GaussianFit3D_Theta.h"

Palm::Palm(const unsigned short* image, double* points, const size_t pointsListSize, const size_t height, const size_t width,
           const size_t waveletNo, const double thresholdVal, const double watershedRatio, const double volMin, const double intMin,
           const size_t gaussFit, const double sigmaX, const double sigmaY, const double theta, const size_t roi)
	: image_(image), points_(points), height_(height), width_(width), size_(static_cast<size_t>(height) * width),
	  nPointsMaxPerPlane_(pointsListSize / N_SEGMENTS), waveletNo_(waveletNo), currentPoint_(0), roi_(roi),
	  thresholdVal_(thresholdVal), volMin_(volMin), intMin_(intMin), watershedRatio_(watershedRatio), gaussFit_(gaussFit),
	  initSigmaX_(sigmaX), initSigmaY_(sigmaY), initTheta_(theta), segmentPub_(height, width)
{
	structFiltParam_.fill(0);
	structFiltParamUse_.fill(0);
	planesToUse_.fill(0);
	operators_.fill(0);
	threshold_.fill(0);

	planesToUse_[waveletNo_] = 2;
	threshold_[waveletNo_] = thresholdVal_;

	if (volMin_ > 0.) {
		structFiltParamUse_[0] = 1; // minimum volume
		structFiltParam_[0] = volMin_;
	}
	if (intMin_ > 0.) {
		structFiltParamUse_[2] = 1; // minimum intensity
		structFiltParam_[2] = intMin_;
	}

	segPub_.resize(size_);
	outputPub_.resize(size_);
	aTrousPub_.wDecomp_.resize(waveletNo + 2);
	for (auto& wavelet : aTrousPub_.wDecomp_) { wavelet.resize(size_); }

	nbPlane_ = 0;
	for (size_t i = 0; i < N_PLANES; ++i) {
		planesToUseB_[i] = planesToUse_[i] != 0;
		if (planesToUseB_[i]) { nbPlane_ = i + 2; }
	}
	planesToUseB_[8] = false;
}

int Palm::Processing()
{
	omp_set_dynamic(1);                                                         // Autorise l’adaptation automatique
	omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max

	int ret;
	if ((ret = aTrousPub_.ATrousDecomposition2D(image_, height_, width_, planesToUseB_.data(), nbPlane_)) < 0) { return ret; }

	for (size_t j = 0; j < size_; ++j) {
		if (aTrousPub_.wDecomp_[waveletNo_][j] > thresholdVal_) {
			outputPub_[j] = aTrousPub_.wDecomp_[waveletNo_][j];
			segPub_[j] = 1;
		}
		else {
			segPub_[j] = 0;
			outputPub_[j] = 0;
		}
	}

	segmentPub_.WatershedIntelligent(outputPub_, segPub_, watershedRatio_);
	segmentPub_.ComputeStatistics(outputPub_, segPub_, aTrousPub_.wDecomp_[waveletNo_], gaussFit_);
	if (volMin_ > 0. || intMin_ > 0.) { segmentPub_.FilterData(segPub_, structFiltParam_.data(), structFiltParamUse_.data()); }

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
				GaussianFit2D fit;
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, offset, rss, roi_);
				sigmaX = initSigmaX_;
				sigmaY = initSigmaX_;
				theta = 0;
				break;
			}

			case 2: { // Gaussian fit of 2D
				GaussianFit2D_Sigma fit;
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, offset, rss, roi_);
				sigmaY = sigmaX;
				theta = 0;
				break;
			}

			case 3: { // Gaussian fit of 3D Theta is constant
				GaussianFit3D fit;
				fit.Run(image_, width_, height_, x, y, intensity, sigmaX, sigmaY, theta, offset, rss, roi_);
				break;
			}

			case 4: { // Gaussian fit of 3D
				GaussianFit3D_Theta fit;
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
