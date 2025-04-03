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
	pointsList_ = points;
	pointsListSize_ = pointsListSize;
	nbSegParams_ = 13; /*nbSegParams_ = 10;*/
	nbPointsMaxPerPlane_ = pointsListSize_ / nbSegParams_;
	height_ = height;
	width_ = width;
	waveletNo_ = waveletNo;
	thresholdVal_ = thresholdVal;
	watershedRatio_ = watershedRatio;
	aTrousPub_ = new ATrous();
	segPub_ = new USrtArray();
	outputPub_ = new FltArray();
	segmentPub_ = new Segment2D();
	volMin_ = volMin;
	intMin_ = intMin;
	gaussFit_ = gaussFit;
	initSigmaX_ = sigmaX;
	initSigmaY_ = sigmaY;
	initTheta_ = theta;
	roi_ = roi;

	for (double& i : structFiltParam_) { i = 0; }
	for (unsigned short& i : structFiltParamUse_) { i = 0; }
	for (size_t& i : planesToUse_) { i = 0; }
	for (size_t& _operator : operators_) { _operator = 0; }
	for (double& i : threshold_) { i = 0; }


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

	input_ = new USrtArray();
	if (input_->alloc(image_, height_, width_) < 0) { return; }
	if (aTrousPub_->wDecomp->alloc(height_, width_, waveletNo + 2) < 0) { return; }
	if (segPub_->realloc(height_, width_) < 0) { return; }
	if (outputPub_->alloc(height_, width_) < 0) { return; }


	planesToUseB_ = new BlArray(9);

	nbPlane_ = 0;
	for (int i = 0; i < 9; i++) {
		(*planesToUseB_)(i) = planesToUse_[i] != 0;
		if ((*planesToUseB_)(i)) { nbPlane_ = i + 2; }
	}
	(*planesToUseB_)(8) = false;
}

long CPALM::Processing()
{
	int ret;
	long nbStruct;
	int nbPointsSign;
	const size_t Size = height_ * width_;
	if ((ret = aTrousPub_->aTrousDecomposition2D(*input_, *aTrousPub_->wDecomp, *planesToUseB_, nbPlane_, false)) < 0) { return ret; }

	for (long j = 0; j < Size; j++) {
		const long offs = waveletNo_ * Size;
		if (aTrousPub_->wDecomp->data[offs + j] > thresholdVal_) {
			outputPub_->data[j] = aTrousPub_->wDecomp->data[offs + j];
			segPub_->data[j] = 1;
		}
		else {
			segPub_->data[j] = 0;
			outputPub_->data[j] = 0;
		}
	}

	if (watershedRatio_ == 1.) { if ((nbStruct = segmentPub_->segmentation(*outputPub_, *segPub_, 1e-8f, true, false)) < 0) { return nbStruct; } }
	else { if ((nbStruct = segmentPub_->watershedIntelligent(*outputPub_, *segPub_, 1e-8f, watershedRatio_)) < 0) { return nbStruct; } }


	if ((nbPointsSign = segmentPub_->computeStatistics(*outputPub_, *segPub_, aTrousPub_->wDecomp->data + waveletNo_ * Size, gaussFit_)) < 0) {
		return nbPointsSign;
	}


	if (volMin_ > 0. || intMin_ > 0.) { if ((nbStruct = segmentPub_->filterData(*segPub_, structFiltParam_, structFiltParamUse_)) < 0) { return ret; } }

	InfoStruct2D** p = (InfoStruct2D**)segmentPub_->getTabInfo();

	for (int i = 1; i <= min(nbStruct, nbPointsMaxPerPlane_ - 2); ++i) {
		switch (gaussFit_) {
			case 0:
			case 5: {
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, 0, 0, 0);
				rss_ = 0;
				fillPointList(p[i]->data);
				break;
			}

			case 1: // Sigma is constant
			{
				// Gaussian fit of 2D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, initSigmaX_, 0, 0);
				GaussianFit2D* fit = new GaussianFit2D();
				fit->Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, offset_, rss_, 10, roi_);
				sigmaY_ = sigmaX_;
				fillPointList(p[i]->data);
				break;
			}

			case 2: {
				// Gaussian fit of 2D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, initSigmaX_, 0, 0);
				GaussianFit2D_Sigma* fit = new GaussianFit2D_Sigma();
				fit->Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, offset_, rss_, 10, roi_);
				sigmaY_ = sigmaX_;
				fillPointList(p[i]->data);
				break;
			}

			case 3: // Theta is constant
			{
				// Gaussian fit of 3D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, initSigmaX_, initSigmaY_, initTheta_);
				GaussianFit3D* fit = new GaussianFit3D();
				fit->Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_, 10, roi_);
				fillPointList(p[i]->data);
				break;
			}

			case 4: {
				// Gaussian fit of 3D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, initSigmaX_, initSigmaY_, initTheta_);
				GaussianFit3D_Theta* fit = new GaussianFit3D_Theta();
				fit->Run(image_, width_, height_, x_, y_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_, 20 * 7, roi_);
				fillPointList(p[i]->data);
				break;
			}

			default: // Error on Gauss Fit Param
				std::fill_n(pointsList_ + currentPoint_, 13, -1);
				currentPoint_ += 13;
				break;
		}
	}

	// Blank Line
	std::fill_n(pointsList_ + currentPoint_, 13, -1);
	currentPoint_ += 13;

	return currentPoint_;
}


CPALM::~CPALM()
{
	delete aTrousPub_;
	delete segPub_;
	delete outputPub_;
	delete segmentPub_;
	delete input_;
	delete planesToUseB_;
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
	pointsList_[currentPoint_++] = sigmaX_;
	pointsList_[currentPoint_++] = sigmaY_;
	pointsList_[currentPoint_++] = theta_;
	pointsList_[currentPoint_++] = y_;
	pointsList_[currentPoint_++] = x_;
	pointsList_[currentPoint_++] = intensity_;
	pointsList_[currentPoint_++] = offset_;
	pointsList_[currentPoint_++] = rss_;
	pointsList_[currentPoint_++] = data.intensity;
	pointsList_[currentPoint_++] = data.surf;

	pointsList_[currentPoint_++] = 0.; // Z coordinate
	pointsList_[currentPoint_++] = 0.; // Spectral or other
	pointsList_[currentPoint_++] = 0.; // Other
}
