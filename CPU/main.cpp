#include "ExportAPI.hpp"
#include "palm/Palm.hpp"
#include "tracking/Parameters.hpp"
#include "tracking/Tracking.hpp"


//-------------------------------------------------------------------------------------------------
EXPORT_API
long CALL_CONV Localization(const unsigned short* image, double* pointsList, const unsigned long pointsListSize, const unsigned long height,
                            const unsigned long width, const unsigned long planes, const double thresholdVal, const double watershedRatio,
                            const unsigned short fit, const double* fit_params)
{
	auto palm = palm_tracer::Palm(image, pointsList, pointsListSize, height, width, planes, thresholdVal, watershedRatio, fit, fit_params);
	return palm.Processing();
}


EXPORT_API
long CALL_CONV Tracking(const double* pointsList, double* trc, const double maxSpeed, const unsigned long minLife,
                        const double decrease, const double costBirth, const unsigned long nbTimes)
{
	const auto Param = new palm_tracer::Parameters();

	Param->maxSpeed_       = maxSpeed;
	Param->costBirthDeath_ = costBirth;

	palm_tracer::Tracking tracking(pointsList, trc, minLife, decrease / 100.0, nbTimes, Param);
	const size_t Count = tracking.Processing();
	delete Param;

	return static_cast<long>(Count);
}
