#include "ExportAPI.hpp"
#include "palm/Palm.hpp"


//-------------------------------------------------------------------------------------------------
EXPORT_API
long CALL_CONV Process(const unsigned short* image, double* pointsList, const unsigned long pointsListSize, const unsigned long height,
                       const unsigned long width, const double thresholdVal, const double watershedRatio, const unsigned short gaussFit,
                       const double sigma, const double angle, const unsigned short roi)
{
	auto palm = Palm(image, pointsList, pointsListSize, height, width, thresholdVal, watershedRatio, gaussFit, sigma, angle, roi);
	return palm.Processing();
}
