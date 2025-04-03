#pragma once
#include <array>

#include "aTrous.h"
#include "segmentation_2D.h"
#include "tools/allDataLib.h"

// Cette classe est exportée de PALM.dll
//class PALM_API CPALM {

constexpr size_t N_SEGMENTS = 13;

class CPALM
{
public:
	CPALM(unsigned short* image, double* points, unsigned long pointsListSize, unsigned long height, unsigned long width,
	      unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX,
	      double sigmaY, double theta, unsigned short roi);
	~CPALM() = default;

	long Processing();

private:
	unsigned short* image_;
	double* points_;
	size_t height_, width_, nPoints_, nPointsMaxPerPlane_, waveletNo_, currentPoint_, roi_;
	double thresholdVal_, volMin_, intMin_, watershedRatio_;

	size_t gaussFit_;
	double initSigmaX_, initSigmaY_, initTheta_;

	std::array<double, N_SEGMENTS> structFiltParam_;
	std::array<unsigned short, N_SEGMENTS> structFiltParamUse_;
	std::array<size_t, 10> planesToUse_;
	std::array<size_t, 9> operators_;
	std::array<double, 9> threshold_;

	ATrous aTrousPub_;
	USrtArray segPub_;
	FltArray outputPub_;
	Segment2D segmentPub_;
	USrtArray input_;
	BlArray planesToUseB_ = BlArray(9);
	size_t nbPlane_;

	// Gaussian Fit Parameters
	double y_, x_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_;

	void initGaussianFitParam(double x, double y, double sigmaX, double sigmaY, double theta);
	void fillPointList(const InfoDataRich2D& data);
};
