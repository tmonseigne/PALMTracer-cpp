#pragma once
#include "aTrous.h"
#include "segmentation_2D.h"
#include "tools/allDataLib.h"

// Cette classe est exportée de PALM.dll
//class PALM_API CPALM {

class CPALM
{
public:
	CPALM(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY,
	      unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX,
	      double sigmaY, double theta, unsigned short roi);
	~CPALM();

	long Processing();

private:

	unsigned long sizeX_;
	unsigned long sizeY_;
	unsigned short* image_;
	double* pointsList_;
	long pointsListSize_;
	int nbPointsMaxPerPlane_;
	int nbSegParams_;
	unsigned long waveletNo_;
	double thresholdVal_;
	double volMin_;
	double intMin_;
	double watershedRatio_;
	unsigned short gaussFit_;
	float initSigmaX_;
	float initSigmaY_;
	float initTheta_;
	unsigned short roi_;
	long currentPoint_;

	double structFiltParam_[13];
	unsigned short structFiltParamUse_[13];
	unsigned short planesToUse_[10];
	unsigned short operators_[9];
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
