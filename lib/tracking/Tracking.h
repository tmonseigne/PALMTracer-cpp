/*******************************************************************************************************************
									Tracking Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/
#pragma once

#include "Parameters.h"
#include "TimeNodes.h"

class Tracking
{
public:
	size_t nTraj_;

	explicit Tracking(const double* points, double* tracks, const size_t minLife, const double decrease, const size_t ntimes, Parameters* params)
		: nTraj_(0), points_(points), tracks_(tracks), decrease_(decrease), minLife_(minLife), nbTimes_(ntimes), params_(params), timeNodes_(nullptr) { }

	virtual ~Tracking();

	TimeNodes* operator()(const size_t i) const { return timeNodes_[i]; }
	TimeNodes** operator[](const size_t i) const { return timeNodes_ + i; }

	void Processing();

protected:
	const double* points_;
	double* tracks_;
	double decrease_;
	size_t minLife_, nbTimes_;
	Parameters* params_;
	TimeNodes** timeNodes_;

	void updateConstantes() const;
	void computeConstantes() const;
	void linkNeighborsBeforeAndAfter() const { for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t]->LinkNeighborsBeforeAndAfter(*timeNodes_[t + 1]); } }
	void linkNeighborsBlinkBeforeAndAfter() const;
	void sortByIntensityAllNeighbors() const { for (size_t t = 0; t < nbTimes_; ++t) { timeNodes_[t]->SortByIntensityAllNeighbors(); } }
	void readPointsList();
	void printConnectionsASCII_trc();
	void setColorToTrackedStructures() const;
	void setIndexToTrackedStructures() const;
	TimeNodes* allocTimeNodes(const size_t time, const size_t nNodes) const { return new TimeNodes(time, nNodes, params_); }
	void allocTabTimes();
};
