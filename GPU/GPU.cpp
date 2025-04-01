#include <cstring>
#include <ctime>
#include <process.h>

#define MAX_THREADS  32

// DLL CPU/GPU for segmentattion
_declspec(dllexport) void _stdcall __OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
_declspec(dllexport) long _stdcall __PALMProcessing();
_declspec(dllexport) void _stdcall __closePALMProcessing();

_declspec(dllexport) void _stdcall InitializeMemory(unsigned short* frames);
_declspec(dllexport) void _stdcall LaunchGaussianFit(double* pointsList, unsigned long pointsListSize, unsigned short sizeStk, unsigned short sizeDataX, unsigned short sizeDataY, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size, unsigned short iterNum, unsigned short fromcplusplus, unsigned short computeZ, double* SigmaXCoefs, double* SigmaYCoefs, double ThreeDThickness/*µm*/, double ThreeDSampling/*µm*/, double pixelSize /*µm*/);


BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
	switch (ul_reason_for_call) {
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH: 
		default: break;
	}
	return TRUE;
}

GPU::GPU(unsigned short* _images, double* _pointsList, unsigned long _pointsListSize, unsigned long _sizeStk, unsigned long _sizeDataX, unsigned long _sizeDataY, unsigned long _waveletNo, double _thresholdVal, double _watershedRatio, double _volMin, double _intMin, unsigned short _gaussFit, double _sigmaX, double _sigmaY, double _angle, unsigned short _size, unsigned short _IsComputingZ, double* _SigmaXCoefs, double* _SigmaYCoefs, double _ThreeDThickness/*µm*/, double _ThreeDSampling/*µm*/, double _pixelSize /*µm*/)
{
	images = _images;
	pointsList = _pointsList;
	pointsListSize = _pointsListSize;
	NbSegParams = 13;
	sizeStk = _sizeStk;
	nbPointsMaxPerPlane = pointsListSize / (NbSegParams * sizeStk);
	sizeDataX = _sizeDataX;
	sizeDataY = _sizeDataY;
	waveletNo = _waveletNo;
	thresholdVal = _thresholdVal;
	watershedRatio = _watershedRatio;
	volMin = _volMin;
	intMin = _intMin;
	gaussFit = _gaussFit;
	sigmaX = _sigmaX;
	sigmaY = _sigmaY;
	angle = _angle;
	sizeFit = _size;

	IsComputingZ = _IsComputingZ;
	SigmaXCoefs = _SigmaXCoefs;
	SigmaYCoefs = _SigmaYCoefs;
	ThreeDThickness = _ThreeDThickness;
	ThreeDSampling = _ThreeDSampling;
	pixelSize = _pixelSize;
}

GPU::~GPU() { }

int GPU::runGPU()
{
	int nb_points = 0;
	unsigned short* image = (unsigned short*)malloc(sizeDataX * sizeDataY * sizeof(unsigned short));
	double* points = (double*)malloc(nbPointsMaxPerPlane * NbSegParams * sizeof(double));

	// 2D detection on CPU/GPU (Live_PALM.dll) version
	long currentPoint = 0;
	for (int plane = 0; plane < sizeStk; plane++) {
		memcpy(image, images + plane * sizeDataX * sizeDataY, sizeDataX * sizeDataY * sizeof(unsigned short));

		__OpenPALMProcessing(image, points, nbPointsMaxPerPlane * NbSegParams, sizeDataY, sizeDataX, waveletNo, thresholdVal, watershedRatio, volMin, intMin, (unsigned short)0, sigmaX, sigmaY, angle, sizeFit);
		int ret_val = __PALMProcessing();
		__closePALMProcessing();

		nb_points += ret_val / NbSegParams - 1;

		// copy data from 'points' to 'pointsList'
		for (int li = 0; li < ret_val / NbSegParams - 1; li++) {
			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = 0.;
			//pointsList[currentPoint++] = points[li * NbSegParams + 3];
			//pointsList[currentPoint++] = points[li * NbSegParams + 4];
			pointsList[currentPoint++] = points[li * NbSegParams + 4];
			pointsList[currentPoint++] = points[li * NbSegParams + 3];
			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = points[li * NbSegParams + 8];
			pointsList[currentPoint++] = points[li * NbSegParams + 9];

			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = 0.;
			pointsList[currentPoint++] = 0.;
		}

		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;
		pointsList[currentPoint++] = -1;

		pointsList[currentPoint++] = -1; // Z coordinate
		pointsList[currentPoint++] = -1; // Spectral or other
		pointsList[currentPoint++] = -1; // Other
	}

	// cleanup memory
	free(image);
	free(points);

	InitializeMemory(images);

	unsigned short iterNum;

	if (gaussFit == 0 || gaussFit == 5) { iterNum = 0; }
	if (gaussFit == 1 || gaussFit == 2 || gaussFit == 3) { iterNum = 10; }
	if (gaussFit == 4) { iterNum = 20 * 7; }

	LaunchGaussianFit(pointsList, pointsListSize, sizeStk, sizeDataY, sizeDataX, gaussFit, sigmaX, sigmaY, angle, sizeFit, iterNum, 1, IsComputingZ, SigmaXCoefs, SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize /*µm*/);

	return currentPoint;
}

GPU* g;

_declspec(dllexport) void _stdcall OpenPALMProcessing(unsigned short* images, double* pointsList, unsigned long pointsListSize, unsigned long sizeStk, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size, unsigned short IsComputingZ, double* SigmaXCoefs, double* SigmaYCoefs, double ThreeDThickness/*µm*/, double ThreeDSampling/*µm*/, double pixelSize /*µm*/) { g = new GPU(images, pointsList, pointsListSize, sizeStk, sizeDataY, sizeDataX, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle, size, IsComputingZ, SigmaXCoefs, SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize /*µm*/); }

_declspec(dllexport) long _stdcall PALMProcessing()
{
	long ret_val = (long)g->runGPU(); // nombre de points détectés  
	return ret_val;
}

_declspec(dllexport) void _stdcall closePALMProcessing() { delete g; }
