#pragma once
#include "palm/ATrous.h"
#include "palm/segmentation_2D.h"
#include "tools/allDataLib.h"

// Cette classe est exportée de PALM.dll
//class PALM_API CPALM {

// includes, project
//#include <helper_cuda.h>
//#include <helper_timer.h>
extern "C"
void cudaTrous2DGPU(float* input, float* output, unsigned int numImages, unsigned short width, unsigned short height, unsigned short step);
extern "C"
void cudaProjectionV02DGPU(unsigned short* input, float* output, unsigned int numImages, unsigned short width, unsigned short height, bool realProj);

extern "C"
void cuda3DFittingGPU_Theta_9x9(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX, float sigmaY, float angle);
extern "C"
void cuda3DFittingGPU_Theta_13x13(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX, float sigmaY, float angle);
extern "C"
void cuda3DFittingGPU_7x7(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX, float sigmaY, float angle, unsigned short computeZ, float* SigmaXCoefs, float* SigmaYCoefs, float ThreeDThickness/*µm*/, float ThreeDSampling/*µm*/, float pixelSize/*µm*/);
extern "C"
void cuda3DFittingGPU_9x9(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX, float sigmaY, float angle, unsigned short computeZ, float* SigmaXCoefs, float* SigmaYCoefs, float ThreeDThickness/*µm*/, float ThreeDSampling/*µm*/, float pixelSize/*µm*/);
extern "C"
void cuda3DFittingGPU_11x11(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX, float sigmaY, float angle, unsigned short computeZ, float* SigmaXCoefs, float* SigmaYCoefs, float ThreeDThickness/*µm*/, float ThreeDSampling/*µm*/, float pixelSize/*µm*/);
extern "C"
void cuda3DFittingGPU_13x13(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX, float sigmaY, float angle, unsigned short computeZ, float* SigmaXCoefs, float* SigmaYCoefs, float ThreeDThickness/*µm*/, float ThreeDSampling/*µm*/, float pixelSize/*µm*/);
extern "C"
void cuda2DFittingGPU_Sigma_7x7(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX);
extern "C"
void cuda2DFittingGPU_Sigma_5x5(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX);
extern "C"
void cuda2DFittingGPU_7x7(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX);
extern "C"
void cuda2DFittingGPU_5x5(unsigned short* images, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int, unsigned short size, float sigmaX);

extern "C"
void GaussianVisu(double* gaussianPointsList, double* pointsList, unsigned short width, unsigned short height, unsigned int nbPoints, unsigned short sizeReg, unsigned short zoom, double conversionRatio /*from: ADU to photons*/, double pixelSize /* unit: nm*/);

extern "C"
void convertUSrtToFloat(unsigned short* input, float* output, unsigned short numImages, unsigned short width, unsigned short height);


class CPALM
{
public:
	// TODO : ajoutez ici vos méthodes.

	long Processing();
	CPALM(unsigned short* image, double* points, unsigned long pointsListSize, unsigned long height, unsigned long width, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double theta, unsigned short roi);
	~CPALM();
	unsigned long height_;
	unsigned long width_;
	unsigned short* image_;
	double* points_;

private:

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
};
