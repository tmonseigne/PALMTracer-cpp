// PALM.cpp : définit le point d'entrée pour l'application DLL.
//


#include "PALM.h"
#include "palm/CPALM.h"

CPALM* c;


_declspec(dllexport)
void _stdcall _OpenPALMProcessing(unsigned short* image, double* pointsList, const unsigned long pointsListSize,
                                  const unsigned long sizeX, const unsigned long sizeY, const unsigned long waveletNo,
                                  const double thresholdVal, const double watershedRatio, const double volMin,
                                  const double intMin, const unsigned short gaussFit, const double sigmaX, const double sigmaY,
                                  const double angle, const unsigned short roi)
{
	c = new CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio,
	              volMin, intMin, gaussFit, sigmaX, sigmaY, angle, roi);
}

_declspec(dllexport)
long _stdcall _PALMProcessing() { return c->Processing(); }

_declspec(dllexport)
void _stdcall _closePALMProcessing() { delete c; }

_declspec(dllexport)
long _stdcall Process(unsigned short* image, double* pointsList, const unsigned long pointsListSize,
                      const unsigned long sizeX, const unsigned long sizeY, const unsigned long waveletNo,
                      const double thresholdVal, const double watershedRatio, const double volMin, const double intMin, const unsigned short gaussFit,
                      const double sigmaX, const double sigmaY, const double angle, const unsigned short roi)
{
	auto palm = CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio,
	                  volMin, intMin, gaussFit, sigmaX, sigmaY, angle, roi);
	return palm.Processing();
}
