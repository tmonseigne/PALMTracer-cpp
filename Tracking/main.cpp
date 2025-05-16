#include "ExportAPI.hpp"
#include "tracking/Parameters.hpp"
#include "tracking/Tracking.hpp"


EXPORT_API
long CALL_CONV Process(const double* pointsList, double* trc, const double maxSpeed, const unsigned long minLife,
                      const double decrease, const double costBirth, const unsigned long nbTimes)
{
	const auto Param = new Parameters();

	Param->maxSpeed_ = maxSpeed;
	Param->costBirthDeath_ = costBirth;

	Tracking tracking(pointsList, trc, minLife, decrease / 100.0, nbTimes, Param);
	tracking.Processing();
	delete Param;

	return static_cast<long>(tracking.nTraj_);
}


EXPORT_API
long CALL_CONV Tracking(const double* points, const unsigned long nPoints, double* trc, const unsigned long nTracks,
                       const double maxSpeed, const double dZdX, const double minLife, const double decrease,
                       const unsigned long allowFusiondisso, const double costBirth, unsigned long nbDim, const unsigned long model,
                       const unsigned long nbTimes)
{
	return Process(points, trc, maxSpeed, static_cast<unsigned long>(minLife), decrease, costBirth, nbTimes);
}
