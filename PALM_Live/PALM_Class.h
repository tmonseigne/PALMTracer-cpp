

#include "allDataLib.h"

#include "aTrous.h"
#include "segmentation_2D.h"

// Cette classe est exportée de PALM.dll
//class PALM_API CPALM {

// includes, project
//#include <helper_cuda.h>
//#include <helper_timer.h>
extern "C"
void cudaTrous2DGPU(float* input, float* output, unsigned int numImages, unsigned short width, unsigned short height, unsigned short step);
extern "C"
void cudaProjectionV02DGPU(unsigned short* input, float* output, unsigned int numImages,unsigned short width, unsigned short height, bool realProj);

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


class CPALM {
public:
	// TODO : ajoutez ici vos méthodes.
	
	long CPALMProcessing();
	CPALM(unsigned short * image,  double * pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, 
							unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin,
							unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
	~CPALM();
	unsigned long _sizeX;
	unsigned long _sizeY;
	unsigned short * _image;
	double * _pointsList;

private:
	unsigned long nbCalls;

	long _pointsListSize;
	int _nbPointsMaxPerPlane;
    int _nbSegParams;
	unsigned long _waveletNo;
	double _thresholdVal;
	double _volMin;
	double _intMin;
	double _watershedRatio;
	unsigned short _gaussFit;
	float _sigmaX;
	float _sigmaY;
	float _angle;
	unsigned short _size;
	long currentPoint;

	double _structFiltParam[13];
	unsigned short _structFiltParamUse[13];
	unsigned short _planesToUse[10];
	unsigned short _operators[9];
	double _threshold[9];

	ATrous *_aTrousPub;
	USrtArray *_segPub;
	FltArray *_outputPub;
	Segment2D *_segmentPub;
	USrtArray *_input;
	BlArray *_planesToUseB;
	BlArray *_planesToUseB2;
	int _nbPlane;
	int waveletSegmentation2D( unsigned short* outputData);
	int waveletSegmentation2D2( unsigned short* outputData);
};