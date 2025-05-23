#include "ExportAPI.hpp"
#include "palm/Palm.hpp"


//-------------------------------------------------------------------------------------------------
EXPORT_API
long CALL_CONV Localisation(const unsigned short* image, double* pointsList, const unsigned long pointsListSize, const unsigned long height,
	const unsigned long width, const unsigned long planes, const double thresholdVal, const double watershedRatio, const unsigned short fit,
                       const double sigma, const double angle, const unsigned short roi)
{
	auto palm = Palm(image, pointsList, pointsListSize, height, width, planes, thresholdVal, watershedRatio, fit, sigma, angle, roi);
	return palm.Processing();
}
