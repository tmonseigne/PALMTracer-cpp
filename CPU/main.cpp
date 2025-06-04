#include <iostream>

#include "ExportAPI.hpp"
#include "alignment/Alignment.hpp"
#include "palm/Palm.hpp"
#include "tracking/Parameters.hpp"
#include "tracking/Tracking.hpp"


//-------------------------------------------------------------------------------------------------
EXPORT_API
long CALL_CONV Localization(const unsigned short* image, double* pointsList, const unsigned long pointsListSize, const unsigned long height,
                            const unsigned long width, const unsigned long planes, const double thresholdVal, const double watershedRatio,
                            const unsigned short fit, const double* fit_params)
{
	const auto Palm = palm_tracer::Palm(image, pointsList, pointsListSize, height, width, planes, thresholdVal, watershedRatio, fit, fit_params);
	return Palm.Processing();
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


EXPORT_API
void CALL_CONV Alignment(const unsigned short* input, unsigned short* output, const unsigned long height,
                         const unsigned long width, const unsigned long planes, const double* factors, const unsigned long upsampling)
{
	std::array<std::array<double, 10>, 2> f;
	size_t idx = 0;
	for (size_t i = 0; i < 2; ++i) { for (size_t j = 0; j < 10; ++j, ++idx) { f[i][j] = factors[idx]; } }

	palm_tracer::Alignment align;
	align.Align(input, output, height, width, planes, f, upsampling);
}
