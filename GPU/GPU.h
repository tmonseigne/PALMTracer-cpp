// includes, project

#include <windows.h>

class GPU
{
public :
	GPU(unsigned short* images, double* pointsList, unsigned long pointsListSize, unsigned long sizeStk, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size, unsigned short IsComputingZ, double* SigmaXCoefs, double* SigmaYCoefs, double ThreeDThickness/*µm*/, double ThreeDSampling/*µm*/, double pixelSize /*µm*/);
	~GPU();
	int runGPU();
	double* pointsList;

private :
	unsigned short* images;

	int nbPointsMaxPerPlane;
	unsigned long pointsListSize;
	unsigned long NbSegParams;
	unsigned long sizeStk;
	unsigned long sizeDataX;
	unsigned long sizeDataY;
	unsigned long waveletNo;
	double thresholdVal;
	double watershedRatio;
	double volMin;
	double intMin;
	unsigned short gaussFit;
	double sigmaX;
	double sigmaY;
	double angle;
	unsigned short sizeFit;

	unsigned short IsComputingZ;
	double ThreeDThickness/*µm*/;
	double ThreeDSampling/*µm*/;
	double* SigmaXCoefs;
	double* SigmaYCoefs;
	double pixelSize;
};
