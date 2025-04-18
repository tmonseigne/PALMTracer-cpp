// Tracking.cpp : Defines the exported functions for the DLL application.
#include <filesystem>
#include "tracking/Parameters.h"
#include "tracking/TabTimes.h"
#include "tracking/DistanceIntensity/TabTimesDistanceIntensity.h"


_declspec(dllexport) long _stdcall Process(double* pointsList, const unsigned long pointsListSize, double* trc, const unsigned long trcSize,
                                           const double maxSpeed, const double dZdX, const unsigned long minLife, const double decrease,
                                           const unsigned long allowFusiondisso, const double costBirth, const unsigned long nbTimes)
{
	const auto Param = new Parameters();

	Param->points_ = pointsList;
	Param->nPoints_ = pointsListSize;
	Param->tracks_ = trc;
	Param->nTracks_ = trcSize;
	Param->maxSpeed_ = maxSpeed;
	Param->dZdx_ = dZdX;
	Param->minLife_ = minLife;
	Param->decrease_ = decrease;
	Param->allowFusionAndFission_ = allowFusiondisso;
	Param->costBirthDeath_ = costBirth;

	Param->nTimes_ = nbTimes;
	Param->nTraj_ = 0;

	TabTimesDistanceIntensity tabTimes(Param);

	int ret;
	if ((ret = tabTimes.Initialisation()) < 0) { return ret; }
	if ((ret = tabTimes.RunModel()) < 0) { return ret; }
	if ((ret = tabTimes.GiveResult()) < 0) { return ret; }

	const long N = static_cast<long>(Param->nTraj_);
	delete Param;

	return N;
}


_declspec(dllexport) long _stdcall Tracking(double* pointsList, const unsigned long pointsListSize, double* trc, const unsigned long trcSize,
                                            const double maxSpeed, const double dZdX, const double minLife, const double decrease,
                                            const unsigned long allowFusiondisso, const double costBirth, unsigned long nbDim, const unsigned long model,
                                            const unsigned long nbTimes)
{
	return Process(pointsList, pointsListSize, trc, trcSize, maxSpeed, dZdX, static_cast<unsigned long>(minLife),
	               decrease, allowFusiondisso, costBirth, nbTimes);
}
