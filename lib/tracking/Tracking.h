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
		: nTraj_(0), points_(points), tracks_(tracks), decrease_(decrease), minLife_(minLife), nbTimes_(ntimes), params_(params) { }

	~Tracking() { timeNodes_.clear(); }

	void Processing();

protected:
	const double* points_;
	double* tracks_;
	double decrease_;
	size_t minLife_, nbTimes_;
	Parameters* params_;
	std::vector<TimeNodes> timeNodes_;

	void updateConstantes() const;
	void computeConstantes() const;
	void linkNeighborsBeforeAndAfter() { for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].LinkNeighborsBeforeAndAfter(timeNodes_[t + 1]); } }
	void linkNeighborsBlinkBeforeAndAfter();
	void sortByIntensityAllNeighbors() { for (size_t t = 0; t < nbTimes_; ++t) { timeNodes_[t].SortByIntensityAllNeighbors(); } }
	void readPointsList();
	void printConnectionsASCII_trc();
	void setColorToTrackedStructures();
	void setIndexToTrackedStructures();
};
