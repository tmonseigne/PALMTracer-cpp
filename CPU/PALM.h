#pragma once

_declspec(dllexport) long _stdcall Process(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY,
                                           unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin,
                                           unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
_declspec(dllexport) void _stdcall _OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX,
                                                       unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin,
                                                       double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
_declspec(dllexport) long _stdcall _PALMProcessing();
_declspec(dllexport) void _stdcall _closePALMProcessing();
