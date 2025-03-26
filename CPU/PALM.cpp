// PALM.cpp : définit le point d'entrée pour l'application DLL.
//


#include "PALM.h"
#include "palm/CPALM.h"

#include "fit/CalibrationFit.h"

#ifdef _MANAGED
#pragma managed(push, off)
#endif

#define min(a,b) (((a) < (b)) ? (a) : (b))

BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
	switch (ul_reason_for_call) {
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH: break;
	}
	return TRUE;
}

#ifdef _MANAGED
#pragma managed(pop)
#endif

CPALM* c;

double sum;
unsigned int pointsNumber;

_declspec(dllexport) void _stdcall _OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
{
	sum = 0.;
	pointsNumber = 0;

	c = new CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle, size);
}

_declspec(dllexport) long _stdcall _PALMProcessing() { return c->CPALMProcessing(); }

_declspec(dllexport) void _stdcall _closePALMProcessing() { delete c; }

_declspec(dllexport) double _stdcall RunCalibrationFit(double* sigma, double* Z, unsigned short pointsNbr, double* C)
{
	float* _sigma = (float*)malloc(pointsNbr * sizeof(float));
	float* _Z = (float*)malloc(pointsNbr * sizeof(float));

	for (int i = 0; i < pointsNbr; i++) {
		_sigma[i] = (float)sigma[i];
		_Z[i] = (float)Z[i];
	}

	float C0 = (float)C[0];
	float C1 = (float)C[1];
	float C2 = (float)C[2];
	float C3 = (float)C[3];
	float C4 = (float)C[4];

	CalibrationFit* fit = new CalibrationFit();
	float RSS = fit->Run(_sigma, _Z, pointsNbr, 150, C0, C1, C2, C3, C4);

	C[0] = (double)C0;
	C[1] = (double)C1;
	C[2] = (double)C2;
	C[3] = (double)C3;
	C[4] = (double)C4;

	// Cleanup memory
	free(_sigma);
	free(_Z);

	return RSS;
}

float GaussianFunction_2D(float x, float y, float I0, float sigma) { return I0 * expf((-powf(x, 2) - powf(y, 2)) / (2 * powf(sigma, 2))); }

_declspec(dllexport) void _stdcall GaussianVisuCPU(double* gaussianPointsList/*(sizeReg*sizeReg*nbPoints)*/, double* pointsList, unsigned long pointsListSize, unsigned short sizeStk, unsigned short width, unsigned short height, unsigned short sizeReg, unsigned short zoom, double conversionRatio, double pixelSize)
{
	unsigned int index = 0;
	float x0, y0;
	unsigned long NbSegParams = 13;

	for (int i = 0; i < (pointsListSize / NbSegParams); i++) {
		if ((int)pointsList[NbSegParams * i + 9] > 0) {
			if ((int)pointsList[NbSegParams * i + 3] != -1) {
				int sizeRegH = sizeReg / 2;

				x0 = (float)pointsList[NbSegParams * i + 4];
				y0 = (float)pointsList[NbSegParams * i + 3];

				x0 = x0 * (float)zoom + 0.5 * ((float)zoom - 1.0);
				y0 = y0 * (float)zoom + 0.5 * ((float)zoom - 1.0);

				float photonsNbr = 2 * 3.14 * (float)pointsList[NbSegParams * i + 5] * (float)pointsList[NbSegParams * i + 0] * (float)pointsList[NbSegParams * i + 0];
				photonsNbr = photonsNbr * conversionRatio;
				float sigma = ((float)pointsList[NbSegParams * i + 0] * pixelSize) / sqrt(photonsNbr); // en nm
				sigma = sigma / (pixelSize / zoom);                                                    // in nm / (pixelSize / zoom); // en pixel SR space
				float A = 10000 / (2 * 3.14 * sigma * sigma);                                          //photonsNbr / (2 * 3.14 * sigma * sigma);

				if (((unsigned int)y0 >= sizeRegH) && ((unsigned int)y0 <= (height * zoom) - 1 - sizeRegH) && ((unsigned int)x0 >= sizeRegH) && ((unsigned int)x0 <= (width * zoom) - 1 - sizeRegH))
					for (int jj = -sizeRegH; jj <= sizeRegH; jj++)
						for (int ii = -sizeRegH; ii <= sizeRegH; ii++) gaussianPointsList[index * sizeReg * sizeReg + (ii + sizeRegH) + (jj + sizeRegH) * sizeReg] = GaussianFunction_2D(x0 - (int)(x0 - ii), y0 - (int)(y0 - jj), A, sigma); // 1000.0f, /*I0*/, 1.0f/*sigma*/);			
			}
			index++;
		}
	}
}
