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
	explicit Tracking(Parameters* params, const size_t ntimes) : nbTimes_(ntimes), parameters_(params), timeNodes_(nullptr) {}
	virtual ~Tracking();
	TimeNodes* operator()(int i) const;
	TimeNodes** operator[](int i) const;


	int Initialisation();
	virtual int RunModel();
	void CheckLinksHomogeneity() const;

protected:
	TimeNodes** timeNodes_;
	void updateConstantes() const;
	void computeConstantes() const;
	int linkNeighborsPresent() const;
	int linkNeighborsBeforeAndAfter() const;
	int linkNeighborsBlinkBeforeAndAfter() const;
	void sortByIntensityAllNeighbors() const;
	int readPointsList();
	void printConnectionsASCII_trc() const;
	void setColorToTrackedStructures() const;
	void setIndexToTrackedStructures() const;
	double getLikenessGlobal() const;
	TimeNodes* allocTimeNodes(const long time, const long nNodes) const { return new TimeNodes(time, nNodes, parameters_); }
	void allocTabTimes();

};
