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
	virtual int RunModel() { return T_ERR_OOP; }
	void CheckLinksHomogeneity() const;

protected:
	TimeNodes** timeNodes_;
	int linkNeighborsPresent() const;
	int linkNeighborsBeforeAndAfter() const;
	int linkNeighborsBlinkBeforeAndAfter() const;
	int sortByIntensityAllNeighbors() const;
	int readPointsList();
	int printConnectionsASCII_trc() const;
	void setColorToTrackedStructures() const;
	void setIndexToTrackedStructures() const;
	virtual TimeNodes* allocTimeNodes(long time, long nbNodes) { return nullptr; }
	void allocTabTimes();

};
