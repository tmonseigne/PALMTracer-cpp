#pragma once

_declspec(dllexport) long _stdcall Process(const unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY,
                                           unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin,
                                           unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short roi);
_declspec(dllexport) void _stdcall _OpenPALMProcessing(const unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX,
                                                       unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin,
                                                       double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short roi);
_declspec(dllexport) long _stdcall _PALMProcessing();
_declspec(dllexport) void _stdcall _closePALMProcessing();
