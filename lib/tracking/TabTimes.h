/*******************************************************************************************************************
									TabTimes Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/
#pragma once

#include "Constants.hpp"
#include "Parameters.h"
#include "TimeNodes.h"

class TabTimes
{
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

public:
	int GiveResult();
	explicit TabTimes(Parameters* publicParameters);
	virtual ~TabTimes();
	TimeNodes* operator()(int i) const;
	TimeNodes** operator[](int i) const;
	int nbTimes_;
	Parameters* parameters_;

	int Initialisation();
	virtual int RunModel() { return T_ERR_OOP; }
	void CheckLinksHomogeneity(const wchar_t* st) const;
};
