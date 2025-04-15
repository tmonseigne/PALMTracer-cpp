#pragma once
#include <array>

#include "ATrous.h"
#include "Segmentation2D.h"

constexpr size_t N_SEGMENTS = 13;
constexpr size_t N_PLANES = 9;

class Palm
{
public:
	Palm(unsigned short* image, double* points, unsigned long pointsListSize, unsigned long height, unsigned long width,
	      unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX,
	      double sigmaY, double theta, unsigned short roi);
	~Palm() = default;

	long Processing();

private:
	unsigned short* image_;
	double* points_;
	size_t height_, width_, size_, nPointsMaxPerPlane_, waveletNo_, currentPoint_, roi_;
	double thresholdVal_, volMin_, intMin_, watershedRatio_;

	size_t gaussFit_;
	double initSigmaX_, initSigmaY_, initTheta_;

	std::array<double, N_SEGMENTS> structFiltParam_;
	std::array<unsigned short, N_SEGMENTS> structFiltParamUse_;
	std::array<size_t, 10> planesToUse_;
	std::array<size_t, N_PLANES> operators_;
	std::array<double, N_PLANES> threshold_;

	ATrous aTrousPub_;
	std::vector<unsigned short> segPub_;
	std::vector<float> outputPub_;
	Segmentation2D segmentPub_;
	std::array<bool, N_PLANES> planesToUseB_;
	size_t nbPlane_;

	// Gaussian Fit Parameters
	double y_, x_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_;

	void initGaussianFitParam(double x, double y, double sigmaX, double sigmaY, double theta);
	void fillPointList(const InfoDataRich2D& data);
};
