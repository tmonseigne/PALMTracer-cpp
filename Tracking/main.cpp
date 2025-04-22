// Tracking.cpp : Defines the exported functions for the DLL application.
#include <filesystem>
#include "tracking/Parameters.h"
#include "tracking/Tracking.h"


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

	Param->nTraj_ = 0;

	Tracking tracking(Param, nbTimes);

	int ret;
	tracking.Initialisation();
	if ((ret = tracking.RunModel()) < 0) { return ret; }
	tracking.GiveResult();

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
