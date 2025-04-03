#pragma once
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
	~CPALM();

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

	double structFiltParam_[13];
	unsigned short structFiltParamUse_[13];
	size_t planesToUse_[10];
	size_t operators_[9];
	double threshold_[9];

	ATrous* aTrousPub_;
	USrtArray* segPub_;
	FltArray* outputPub_;
	Segment2D* segmentPub_;
	USrtArray* input_;
	BlArray* planesToUseB_;
	int nbPlane_;

	// Gaussian Fit Parameters
	double y_, x_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_;

	void initGaussianFitParam(double x, double y, double sigmaX, double sigmaY, double theta);
	void fillPointList(const InfoDataRich2D& data);
};
