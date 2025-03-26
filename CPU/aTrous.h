#ifndef __class_a_trous__
#define __class_a_trous__

//#include <iostream.h>
#include <memory.h>
#include "tools/allDataLib.h"
#include "tools/miscTools.h"


class ATrous
{
	//general part:
private:
	FltArray* temp1;
	FltArray* temp2;
	FltArray* inputV0;
	FltArray* VN;
	FltArray* VNplus1;

public:
	FltArray* wDecomp;
	ATrous();
	~ATrous();

	//tools:
public:
	double stdRef(long dim, unsigned long realProj, unsigned long plan);

	//1D transform:
private:
	int offsetTest(int i, int maxSize);
	int testSizesWD1D(FltArray& waveletDecomp, FltArray& reconstruction, BlArray& planesToUse, unsigned long& nbPlane);
	int testSizes1D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	int testSizes1D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	void projectionV01D(unsigned short* srcData, float* destData, unsigned long realProj, unsigned long sizeData);
	void projectionV01D(float* srcData, float* destData, unsigned long realProj, unsigned long sizeData);
	int b3SplineFiltering1D(float* srcData, float* destData, unsigned long sizeData, unsigned long currentPlane);

public:
	int aTrousFiltering1D(FltArray& inputData, FltArray& outputData, unsigned short* planesToUse, double* threshold, unsigned long nbPlane, bool realProj, bool useLastPlane, bool adjoint);
	int aTrousDecomposition1D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
	int aTrousDecomposition1D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
	int aTrousReconstruction1D(FltArray& waveletDecomp, FltArray& reconstruction, BlArray& planesToUse, unsigned long nbPlane, bool useLastPlane, bool adjoint);

	//2D transform:
private:
	int offsetTest(int i, int j, int maxSizeI, int maxSizeJ);
	int offset(int i, int j, int maxSizeJ);
	//int testSizesWD2D(FltArray &waveletDecomp, FltArray &reconstruction, BlArray &planesToUse, unsigned long &nbPlane);
	int testSizes2D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	int testSizes2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane);
	void projectionV02D(unsigned short* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY);
	void projectionV02D(float* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY);
	int b3SplineFiltering2D(float* srcData, float* destData, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long currentPlane);
	void threshold(FltArray& WAngle, double thre, long currentPlane);
	//int directional2B3SplineFiltering2D(FltArray &_WNplus1, FltArray &_VN,unsigned long currentPlane, long indexAngle);
	int directional3B3SplineFiltering2D(FltArray& _WNplus1, FltArray& _VN, unsigned long currentPlane, long indexAngle);
	int directional3B3SplineFiltering2D(FltArray& _WNplus1, FltArray& _VN, unsigned long currentPlane, long indexAngle, USrtArray& mask, bool unsample);
	//int directionalB3SplineFilteringTotal2D(FltArray &_WNplus1, FltArray &_VNplus1, FltArray &_VN,unsigned long currentPlane, double angle, 
	//										 unsigned long indiceAngle, unsigned long nbAngle, unsigned long transformType);
	int directionalB3SplineFiltering2D(FltArray& _WNplus1, FltArray& _VNplus1, FltArray& _VN, unsigned long currentPlane, double angle, unsigned long transformType);

public:
	int aTrousDecomposition2D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
	int aTrousDecomposition2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj);
	int aTrousReconstruction2D(FltArray& waveletDecomp, FltArray& reconstruction, BlArray& planesToUse, unsigned long nbPlane, bool useLastPlane, bool adjoint);
	int aTrousDecompAniso2D(USrtArray& inputData, FltArray& waveletPlanes, BlArray& planesToUseB, unsigned short* planesToUse, double* threshold, unsigned long nbPlane, unsigned long nbAngle, bool realProj, bool useLastPlane, unsigned long transformType, short autoThresholdModel);
	//int aTrousDecompAnisoTotal2D( USrtArray &inputData, FltArray &waveletPlanes, BlArray &planesToUseB, 
	//					unsigned short* planesToUse, double *threshold, unsigned long nbPlane, 
	//					unsigned long nbAngle, bool realProj, bool useLastPlane, unsigned long transformType);
	int aTrousDecompAniso2DHardCTM(USrtArray& inputData, FltArray& waveletPlanes, FltArray& filteringData, BlArray& planesToUseB, unsigned short* planesToUse, double* threshold, unsigned long nbPlane, unsigned long nbAngle, bool realProj, bool useLastPlane, unsigned long transformType);
	int aTrousEnergiesAniso2D(FltArray& VnMap, FltArray& directionnalVnp1, USrtArray& mask, unsigned long currentPlane, bool unsample);

	//3D transform
private:
	int offsetTest(int i, int j, int k, int maxSizeI, int maxSizeJ, int maxSizeK);
	int offset(int i, int j, int k, int maxSizeJ, int maxSizeK);
	int testSizesWD3D(FltArray& waveletDecomp, FltArray& reconstruction, BlArray& planesToUse, unsigned long& nbPlane);
	int testSizes3D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane, double deltaPlan);
	int testSizes3D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane, double deltaPlan);
	void projectionV03D(unsigned short* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long sizeDataZ);
	void projectionV03D(float* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long sizeDataZ);
	int b3SplineFiltering3D(float* srcData, float* destData, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long sizeDataZ, unsigned long currentPlane, double deltaPlan);
	int b3SplineFiltering3D(float* srcData, float* destData, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long sizeDataZ, unsigned long currentPlane);

public:
	int aTrousDecomposition3D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj, double deltaPlan);
	int aTrousDecomposition3D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj, double deltaPlan);
	int aTrousReconstruction3D(FltArray& waveletDecomp, FltArray& reconstruction, BlArray& planesToUse, unsigned long nbPlane, bool useLastPlane, bool adjoint, double deltaPlan);
};


#endif //__class_a_trous__
