#include "ExportAPI.hpp"
#include "palm/Palm.hpp"

Palm* c;

EXPORT_API
long CALL_CONV Process(const uint16_t* image, double* pointsList, const uint32_t pointsListSize,
                       const uint32_t sizeX, const uint32_t sizeY, const uint32_t waveletNo,
                       const double thresholdVal, const double watershedRatio, const double volMin, const double intMin, const uint16_t gaussFit,
                       const double sigmaX, const double sigmaY, const double angle, const uint16_t roi)
{
	auto palm = Palm(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio,
	                 volMin, intMin, gaussFit, sigmaX, sigmaY, angle, roi);
	return palm.Processing();
}
