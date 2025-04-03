// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "CPALM.h"
#include "fit/GaussianFit2D.h"
#include "fit/GaussianFit2D_Sigma.h"
#include "fit/GaussianFit3D.h"
#include "fit/GaussianFit3D_Theta.h"

CPALM::CPALM(unsigned short* image, double* points, const unsigned long pointsListSize, const unsigned long height, const unsigned long width,
             const unsigned long waveletNo, const double thresholdVal, const double watershedRatio, const double volMin, const double intMin,
             const unsigned short gaussFit, const double sigmaX, const double sigmaY, const double theta, const unsigned short roi)
{
	currentPoint_ = 0;
	image_ = image;
	points_ = points;
	nPoints_ = pointsListSize;
	nPointsMaxPerPlane_ = nPoints_ / N_SEGMENTS;
	height_ = height;
	width_ = width;
	waveletNo_ = waveletNo;
	thresholdVal_ = thresholdVal;
	watershedRatio_ = watershedRatio;
	volMin_ = volMin;
	intMin_ = intMin;
	gaussFit_ = gaussFit;
	initSigmaX_ = sigmaX;
	initSigmaY_ = sigmaY;
	initTheta_ = theta;
	roi_ = roi;

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

	if (input_.alloc(image_, height_, width_) < 0) { return; }
	if (aTrousPub_.wDecomp_.alloc(height_, width_, waveletNo + 2) < 0) { return; }
	if (segPub_.realloc(height_, width_) < 0) { return; }
	if (outputPub_.alloc(height_, width_) < 0) { return; }


	nbPlane_ = 0;
	for (size_t i = 0; i < N_PLANES; i++) {
		planesToUseB_[i] = planesToUse_[i] != 0;
		if (planesToUseB_[i]) { nbPlane_ = i + 2; }
	}
	planesToUseB_[8] = false;
}

long CPALM::Processing()
{
	int ret;
	long nbStruct;
	int nbPointsSign;
	const size_t Size = height_ * width_;
	if ((ret = aTrousPub_.ATrousDecomposition2D(input_, aTrousPub_.wDecomp_, planesToUseB_.data(), nbPlane_, false)) < 0) { return ret; }

	for (long j = 0; j < Size; j++) {
		const long offs = waveletNo_ * Size;
		if (aTrousPub_.wDecomp_.data[offs + j] > thresholdVal_) {
			outputPub_.data[j] = aTrousPub_.wDecomp_.data[offs + j];
			segPub_.data[j] = 1;
		}
		else {
			segPub_.data[j] = 0;
			outputPub_.data[j] = 0;
		}
	}

	if ((nbStruct = segmentPub_.watershedIntelligent(outputPub_, segPub_, 1e-8f, watershedRatio_)) < 0) { return nbStruct; }


	if ((nbPointsSign = segmentPub_.computeStatistics(outputPub_, segPub_, aTrousPub_.wDecomp_.data + waveletNo_ * Size, gaussFit_)) < 0) {
		return nbPointsSign;
	}


	if (volMin_ > 0. || intMin_ > 0.) { if ((nbStruct = segmentPub_.filterData(segPub_, structFiltParam_.data(), structFiltParamUse_.data())) < 0) { return ret; } }

	auto p = reinterpret_cast<InfoStruct2D**>(segmentPub_.getTabInfo());

	for (int i = 1; i <= min(nbStruct, nPointsMaxPerPlane_ - 2); ++i) {
		initGaussianFitParam(p[i]->data_.moment.y, p[i]->data_.moment.x, initSigmaX_, initSigmaY_, initTheta_);
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

			case 3: {// Theta is constant
			
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
		fillPointList(p[i]->data_);
	}

	// Blank Line
	std::fill_n(points_ + currentPoint_, N_SEGMENTS, -1);
	currentPoint_ += N_SEGMENTS;

	return static_cast<long>(currentPoint_);
}


void CPALM::initGaussianFitParam(const double x, const double y, const double sigmaX, const double sigmaY, const double theta)
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

void CPALM::fillPointList(const InfoDataRich2D& data)
{
	points_[currentPoint_++] = sigmaX_;
	points_[currentPoint_++] = sigmaY_;
	points_[currentPoint_++] = theta_;
	points_[currentPoint_++] = y_;
	points_[currentPoint_++] = x_;
	points_[currentPoint_++] = intensity_;
	points_[currentPoint_++] = offset_;
	points_[currentPoint_++] = rss_;
	points_[currentPoint_++] = data.intensity;
	points_[currentPoint_++] = data.surf;
	points_[currentPoint_++] = 0.; // Z coordinate
	points_[currentPoint_++] = 0.; // Spectral or other
	points_[currentPoint_++] = 0.; // Other
}
