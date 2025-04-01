// PALM.cpp : définit le point d'entrée pour l'application DLL.
//


#include "PALM.h"
#include "palm/CPALM.h"

CPALM* c;


_declspec(dllexport) void _stdcall _OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX,
                                                       unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin,
                                                       double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
{
	c = new CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle,
	              size);
}

_declspec(dllexport) long _stdcall _PALMProcessing() { return c->CPALMProcessing(); }

_declspec(dllexport) void _stdcall _closePALMProcessing() { delete c; }

_declspec(dllexport) long _stdcall Process(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX,
                                           unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin,
                                           double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
{
	auto palm = CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle,
	                  size);
	return palm.CPALMProcessing();
}
