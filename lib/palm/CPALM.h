#pragma once
#include <array>

#include "aTrous.h"
#include "segmentation_2D.h"
#include "tools/allDataLib.h"

// Cette classe est exportée de PALM.dll
//class PALM_API CPALM {

class CPALM
{
public:
	CPALM(unsigned short* image, double* points, unsigned long pointsListSize, unsigned long height, unsigned long width,
	      unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX,
	      double sigmaY, double theta, unsigned short roi);
	~CPALM() = default;

	long Processing();

private:
	size_t height_;
	size_t width_;
	unsigned short* image_;
	double* pointsList_;
	long pointsListSize_;
	int nbPointsMaxPerPlane_;
	int nbSegParams_;
	size_t waveletNo_;
	double thresholdVal_;
	double volMin_;
	double intMin_;
	double watershedRatio_;
	size_t gaussFit_;
	double initSigmaX_;
	double initSigmaY_;
	double initTheta_;
	size_t roi_;
	long currentPoint_;

	std::array<double, 13> structFiltParam_;
	std::array<unsigned short, 13> structFiltParamUse_;
	std::array<size_t, 10> planesToUse_;
	std::array<size_t, 9> operators_;
	std::array<double, 9> threshold_;

	ATrous aTrousPub_;
	USrtArray segPub_;
	FltArray outputPub_;
	Segment2D segmentPub_;
	USrtArray input_;
	BlArray planesToUseB_ = BlArray(9);
	int nbPlane_;

	// Gaussian Fit Parameters
	double y_, x_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_;

	void initGaussianFitParam(double x, double y, double sigmaX, double sigmaY, double theta);
	void fillPointList(const InfoDataRich2D& data);
};
