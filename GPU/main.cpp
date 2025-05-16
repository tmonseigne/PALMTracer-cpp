#include "ExportAPI.hpp"
#include "palm/Palm.hpp"

Palm* c;

EXPORT_API
long CALL_CONV Process(const unsigned short* image, double* pointsList, const unsigned long pointsListSize,
                      const unsigned long sizeX, const unsigned long sizeY, const unsigned long waveletNo,
                      const double thresholdVal, const double watershedRatio, const double volMin, const double intMin, const unsigned short gaussFit,
                      const double sigmaX, const double sigmaY, const double angle, const unsigned short roi)
{
	auto palm = Palm(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio,
	                 volMin, intMin, gaussFit, sigmaX, sigmaY, angle, roi);
	return palm.Processing();
}
