///-------------------------------------------------------------------------------------------------
/// \file TrackAnalysis.cpp
/// <summary>Orchestration et sérialisation des analyses demandées par la DLL.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#include "TrackAnalysis.hpp"

#include "signal_processing/Matrix.hpp"
#include "tracking/Analysis.hpp"
#include "tracking/Tracks.hpp"

namespace palm_tracer::cpu {
bool AnalyzeTracks(const double* input, double* outputMSD, double* outputInd, double* outputFit, const uint64_t nRow,
                   const bool isMSD, const bool isInd, const bool is3D,
                   const double pixelSize, const double exposureTime, const uint64_t fitMode, const double* fitParams)
{
	palm_tracer::tracking::ComputesParams p = {
		.mode = palm_tracer::tracking::fit_mode::linear,
		.fitLength = static_cast<uint64_t>(fitParams[0]),
		.is3D = is3D,
		.pixelSize = pixelSize,
		.exposureTime = exposureTime
	};
	palm_tracer::tracking::Tracks tracks;
	tracks.LoadTracks(input, nRow, pixelSize);

	if (isMSD) {
		const auto MSD = palm_tracer::tracking::analysis::ComputeAllMSD(tracks.GetTracks());
		palm_tracer::matrix::SendIrregularArray(MSD, outputMSD);
	}
	if (isInd) {
		const auto Ind = palm_tracer::tracking::analysis::ComputeAllInstantDiffusion(tracks.GetTracks(), p);
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
		const auto Fit = palm_tracer::tracking::analysis::ComputeAllFitMetrics(tracks.GetTracks(), p);
		palm_tracer::matrix::SendIrregularArray(Fit, outputFit);
	}
	return true;
}
} // namespace palm_tracer::cpu
