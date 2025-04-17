// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "Palm.h"
#include "fit/GaussianFit2D.h"
#include "fit/GaussianFit2D_Sigma.h"
#include "fit/GaussianFit3D.h"
#include "fit/GaussianFit3D_Theta.h"

Palm::Palm(unsigned short* image, double* points, const size_t pointsListSize, const size_t height, const size_t width,
           const size_t waveletNo, const double thresholdVal, const double watershedRatio, const double volMin, const double intMin,
           const size_t gaussFit, const double sigmaX, const double sigmaY, const double theta, const size_t roi)
	: image_(image), points_(points), height_(height), width_(width), size_(static_cast<size_t>(height) * width),
	  nPointsMaxPerPlane_(pointsListSize / N_SEGMENTS), waveletNo_(waveletNo), currentPoint_(0), roi_(roi), thresholdVal_(thresholdVal),
	  volMin_(volMin), intMin_(intMin), watershedRatio_(watershedRatio), gaussFit_(gaussFit), initSigmaX_(sigmaX), initSigmaY_(sigmaY),
	  initTheta_(theta), segmentPub_(height, width), y_(0), x_(0), intensity_(0), sigmaX_(0), sigmaY_(0), theta_(0), offset_(0), rss_(0)
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
	int ret;
	if ((ret = aTrousPub_.ATrousDecomposition2D(image_, height_, width_, planesToUseB_.data(), nbPlane_, false)) < 0) { return ret; }

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

	segmentPub_.WatershedIntelligent(outputPub_, segPub_, 1e-8, watershedRatio_);
	segmentPub_.ComputeStatistics(outputPub_, segPub_, aTrousPub_.wDecomp_[waveletNo_], gaussFit_);
	if (volMin_ > 0. || intMin_ > 0.) { segmentPub_.FilterData(segPub_, structFiltParam_.data(), structFiltParamUse_.data()); }

	const auto& P = segmentPub_.GetInfos();

	for (size_t i = 0; i < std::min<size_t>(P.size(), nPointsMaxPerPlane_ - 1); ++i) {
		initGaussianFitParam(P[i].moment_.y, P[i].moment_.x, initSigmaX_, initSigmaY_, initTheta_);
		switch (gaussFit_) {
			default: {
				sigmaX_ = 0;
				sigmaY_ = 0;
				theta_ = 0;
				rss_ = 0;
				break;
			}

			case 1: { // Sigma is constant

				// Gaussian fit of 2D
				GaussianFit2D fit;
				fit.Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, offset_, rss_, 10, roi_);
				sigmaX_ = initSigmaX_;
				sigmaY_ = initSigmaX_;
				theta_ = 0;
				break;
			}

			case 2: {
				// Gaussian fit of 2D
				GaussianFit2D_Sigma fit;
				fit.Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, offset_, rss_, 10, roi_);
				sigmaY_ = sigmaX_;
				theta_ = 0;
				break;
			}

			case 3: { // Theta is constant

				// Gaussian fit of 3D
				GaussianFit3D fit;
				fit.Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_, 10, roi_);
				break;
			}

			case 4: {
				// Gaussian fit of 3D
				GaussianFit3D_Theta fit;
				fit.Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_, 140, roi_);
				break;
			}
		}
		fillPointList(P[i]);
	}

	// Blank Line
	std::fill_n(points_ + currentPoint_, N_SEGMENTS, -1);
	currentPoint_ += N_SEGMENTS;

	return static_cast<int>(currentPoint_);
}


void Palm::initGaussianFitParam(const double x, const double y, const double sigmaX, const double sigmaY, const double theta)
{
	x_ = x;
	y_ = y;
	intensity_ = 0;
	sigmaX_ = sigmaX;
	sigmaY_ = sigmaY;
	theta_ = theta;
	offset_ = 0;
	rss_ = -1;
}

void Palm::fillPointList(const Info2D& data)
{
	points_[currentPoint_++] = sigmaX_;
	points_[currentPoint_++] = sigmaY_;
	points_[currentPoint_++] = theta_;
	points_[currentPoint_++] = y_;
	points_[currentPoint_++] = x_;
	points_[currentPoint_++] = intensity_;
	points_[currentPoint_++] = offset_;
	points_[currentPoint_++] = rss_;
	points_[currentPoint_++] = data.intensity_;
	points_[currentPoint_++] = data.surf_;
	points_[currentPoint_++] = 0.; // Z coordinate
	points_[currentPoint_++] = 0.; // Spectral or other
	points_[currentPoint_++] = 0.; // Other
}
