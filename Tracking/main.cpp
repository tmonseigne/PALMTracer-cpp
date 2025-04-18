// Tracking.cpp : Defines the exported functions for the DLL application.
#include <filesystem>
#include <string>
#include "Constants.hpp"
#include "tracking/Parameters.h"
#include "tracking/tab_times.h"
#include "tracking/closerNeighbors/tab_times_closer_neighbors.h"
#include "tracking/likenessDistInten/tab_times_dist_inten.h"


_declspec(dllexport) long _stdcall Process(double* pointsList, const unsigned long pointsListSize, double* trc, const unsigned long trcSize,
                                           const double maxSpeed, const double dZdX, const unsigned long minLife, const double decrease,
                                           const unsigned long allowFusiondisso, const double costBirth, const unsigned long model,
                                           const unsigned long nbTimes)
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
	Param->model_ = model;

	Param->nTimes_ = nbTimes;
	Param->nTraj_ = 0;

	TabTimes* tabTimes;
	switch (Param->model_) {
		case 1: tabTimes = new TabTimesCloserNeighbors(Param);
			break;
		case 2: tabTimes = new TabTimesDistanceIntensity(Param);
			break;
		default: return T_ERR_BAD_MODEL;
	}

	int ret;
	if ((ret = tabTimes->initialisation()) < 0) { return ret; }
	if ((ret = tabTimes->runModel()) < 0) { return ret; }
	if ((ret = tabTimes->giveResult()) < 0) { return ret; }

	const long N = static_cast<long>(Param->nTraj_);
	delete Param;
	delete tabTimes;

	return N;
}


_declspec(dllexport) long _stdcall Tracking(double* pointsList, const unsigned long pointsListSize, double* trc, const unsigned long trcSize,
                                            const double maxSpeed, const double dZdX, const double minLife, const double decrease,
                                            const unsigned long allowFusiondisso, const double costBirth, unsigned long nbDim, const unsigned long model,
                                            const unsigned long nbTimes)
{
	return Process(pointsList, pointsListSize, trc, trcSize, maxSpeed, dZdX, minLife, decrease, allowFusiondisso, costBirth, model, nbTimes);
}
