/*******************************************************************************************************************
									Tracking Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/
#pragma once

#include "Constants.hpp"
#include "Parameters.h"
#include "TimeNodes.h"

class Tracking
{
public:
	int nbTimes_;
	Parameters* parameters_;

	int GiveResult();
	explicit Tracking(Parameters* params);
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
	int sortByIntensityAllNeighbors() const;
	int readPointsList();
	int printConnectionsASCII_trc() const;
	void setColorToTrackedStructures() const;
	void setIndexToTrackedStructures() const;
	double getLikenessGlobal() const;
	TimeNodes* allocTimeNodes(const long time, const long nNodes) { return new TimeNodes(time, nNodes, parameters_); }
	void allocTabTimes();

};
