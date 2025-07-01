#include <cstdint>
#include <iostream>

#include "ExportAPI.hpp"
#include "alignment/Alignment.hpp"
#include "palm/Palm.hpp"
#include "tracking/Parameters.hpp"
#include "tracking/Tracking.hpp"


//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
void CALL_CONV Alignment(const uint16_t* input, uint16_t* output, const uint32_t height,
                         const uint32_t width, const uint32_t planes, const double* factors, const uint32_t upsampling)
{
	std::array<std::array<double, 10>, 2> f{};
	size_t idx = 0;
	for (size_t i = 0; i < 2; ++i) { for (size_t j = 0; j < 10; ++j, ++idx) { f[i][j] = factors[idx]; } }

	palm_tracer::Alignment align;
	align.Align(input, output, height, width, planes, f, upsampling);
}

//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
double CALL_CONV AutoThreshold(const uint16_t* image, const uint32_t height, const uint32_t width, const double* fit_params)
{
	const auto PointsListSize = static_cast<size_t>(0.2 * static_cast<double>(height * width) * palm_tracer::N_LOCALIZATIONS);
	std::vector<double> pointsList(PointsListSize);
	auto palm = palm_tracer::Palm(image, pointsList.data(), PointsListSize, height, width, 1, 0, 10, 0, fit_params);
	return palm.AutoThreshold();
}

//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
uint32_t CALL_CONV Localization(const uint16_t* image, double* pointsList, const uint32_t pointsListSize, const uint32_t height,
                                const uint32_t width, const uint32_t planes, const double thresholdVal, const double watershedRatio,
                                const uint16_t fit, const double* fit_params)
{
	const auto Palm = palm_tracer::Palm(image, pointsList, pointsListSize, height, width, planes, thresholdVal, watershedRatio, fit, fit_params);
	return Palm.Processing();
}


//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
uint32_t CALL_CONV Tracking(const double* pointsList, double* trc, const double maxSpeed, const uint32_t minLife,
                            const double decrease, const double costBirth, const uint32_t nbTimes)
{
	const auto Param = new palm_tracer::Parameters();

	Param->maxSpeed_       = maxSpeed;
	Param->costBirthDeath_ = costBirth;

	palm_tracer::Tracking tracking(pointsList, trc, minLife, decrease / 100.0, nbTimes, Param);
	const size_t Count = tracking.Processing();
	delete Param;

	return static_cast<uint32_t>(Count);
}
