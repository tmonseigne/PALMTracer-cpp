#include <cstdint>
#include <iostream>

#include "ExportAPI.hpp"
#include "Palm.hpp"
#include "signal_processing/Alignment.hpp"
#include "signal_processing/Matrix.hpp"
#include "tracking/Parameters.hpp"
#include "tracking/Tracking.hpp"
#include "tracking/Tracks.hpp"


//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
void CALL_CONV Alignment(const uint16_t* input, uint16_t* output, const uint32_t height,
                         const uint32_t width, const uint32_t planes, const double* factors, const uint32_t upsampling)
{
	std::array<std::array<double, 10>, 2> f{};
	size_t idx = 0;
	for (size_t i = 0; i < 2; ++i) { for (size_t j = 0; j < 10; ++j, ++idx) { f[i][j] = factors[idx]; } }

	palm_tracer::Alignment align;
	const palm_tracer::Vec3I Dim = {.x = width, .y = height, .z = planes};
	align.Align(input, output, Dim, f, upsampling);
}

//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
double CALL_CONV AutoThreshold(const uint16_t* image, const uint32_t height, const uint32_t width, const double* fit_params)
{
	const auto PointsListSize = static_cast<size_t>(0.2 * static_cast<double>(height * width) * palm_tracer::N_LOCALIZATIONS);
	std::vector<double> pointsList(PointsListSize);
	const palm_tracer::Vec3I Dim = {.x = width, .y = height, .z = 1};

	auto palm = palm_tracer::Palm(image, pointsList.data(), PointsListSize, Dim, 0, 10, 0, fit_params);
	return palm.AutoThreshold();
}

//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
uint32_t CALL_CONV Localization(const uint16_t* image, double* pointsList, const uint32_t pointsListSize, const uint32_t height,
                                const uint32_t width, const uint32_t planes, const double thresholdVal, const double watershedRatio,
                                const uint16_t fit, const double* fit_params)
{
	const palm_tracer::Vec3I Dim = {.x = width, .y = height, .z = planes};

	const auto Palm = palm_tracer::Palm(image, pointsList, pointsListSize, Dim, thresholdVal, watershedRatio, fit, fit_params);
	return Palm.Processing();
}


//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
uint32_t CALL_CONV Tracking(const double* pointsList, double* trc, const double maxSpeed, const uint32_t minLife,
                            const double decrease, const double costBirth, const uint32_t nbTimes)
{
	const auto Param = new palm_tracer::tracking::Parameters();

	Param->maxSpeed_       = maxSpeed;
	Param->costBirthDeath_ = costBirth;

	palm_tracer::tracking::Tracking tracking(pointsList, trc, minLife, decrease / 100.0, nbTimes, Param);
	const size_t Count = tracking.Processing();
	delete Param;

	return static_cast<uint32_t>(Count);
}

//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
uint32_t CALL_CONV BlinkingReconnection(const double* input, double* output, const uint32_t nRow, const double pixelSize,
                                        const uint32_t mode, const uint32_t maxDuration, const double maxSpeed)
{
	palm_tracer::tracking::Tracks tracks;
	tracks.LoadTracking(input, nRow, pixelSize);
	tracks.BlinkingReconnection(mode, maxDuration, maxSpeed);
	return static_cast<uint32_t>(tracks.FillOutput(output));
}

//-------------------------------------------------------------------------------------------------
extern "C" EXPORT_API
bool CALL_CONV TracksCompute(const double* input, double* outputMSD, double* outputInd, double* outputFit, const uint32_t nRow,
                             const bool isMSD, const bool isInd, const bool is3D,
                             const double pixelSize, const double exposureTime, const uint32_t fitMode, const double* fitParams)
{
	palm_tracer::tracking::ComputesParams p = {
		.mode = palm_tracer::tracking::fit_mode::linear,
		.fitLength = static_cast<size_t>(fitParams[0]),
		.is3D = is3D,
		.pixelSize = pixelSize,
		.exposureTime = exposureTime
	};
	palm_tracer::tracking::Tracks tracks;
	tracks.LoadTracking(input, nRow, pixelSize);

	if (isMSD) {
		const auto MSD = tracks.ComputeAllMSD();
		palm_tracer::matrix::SendIrregularArray(MSD, outputMSD);
	}
	if (isInd) {
		const auto Ind = tracks.ComputeAllInstantDiffusion(p);
		palm_tracer::matrix::SendIrregularArray(Ind, outputInd);
	}
	if (fitMode != 0) {
		switch (fitMode) {
			case 2: {
				p.mode = palm_tracer::tracking::fit_mode::power;
				break;
			}
			case 3: {
				p.mode = palm_tracer::tracking::fit_mode::exponential;
				break;
			}
			default: break;
		}
		const auto Fit = tracks.ComputeAllFitMetrics(p);
		palm_tracer::matrix::SendIrregularArray(Fit, outputFit);
	}
	return true;
}
