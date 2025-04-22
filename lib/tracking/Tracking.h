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
	size_t nbTimes_;
	Parameters* parameters_;

	void GiveResult() const;
	explicit Tracking(Parameters* params, const size_t ntimes) : nbTimes_(ntimes), parameters_(params), timeNodes_(nullptr) { }
	virtual ~Tracking();
	TimeNodes* operator()(const size_t i) const { return timeNodes_[i]; }
	TimeNodes** operator[](const size_t i) const { return timeNodes_ + i; }

	void Initialisation();
	virtual int RunModel();

protected:
	TimeNodes** timeNodes_;
	void updateConstantes() const;
	void computeConstantes() const;
	void linkNeighborsPresent() const { for (size_t t = 0; t < nbTimes_; ++t) { timeNodes_[t]->LinkNeighborsPresent(); } }
	void linkNeighborsBeforeAndAfter() const { for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t]->LinkNeighborsBeforeAndAfter(*timeNodes_[t + 1]); } }
	void linkNeighborsBlinkBeforeAndAfter() const;
	void sortByIntensityAllNeighbors() const;
	void readPointsList();
	void printConnectionsASCII_trc() const;
	void setColorToTrackedStructures() const;
	void setIndexToTrackedStructures() const;
	double getLikenessGlobal() const;
	TimeNodes* allocTimeNodes(const size_t time, const size_t nNodes) const { return new TimeNodes(time, nNodes, parameters_); }
	void allocTabTimes();
};
