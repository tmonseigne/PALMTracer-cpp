/*******************************************************************************************************************
									TabTimes Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/
#pragma once

#include <stdlib.h>
#include <vector>
#include "general.h"
#include "time_nodes.h"
#include "tools/miscTools.h"
//#include "tracking_2D.h"
using namespace std;

class TabTimes
{
protected:
	TimeNodes** timeNodes;
	int linkNeighborsPresent();
	int linkNeighborsBeforeAndAfter();
	int linkNeighborsBlinkBeforeAndAfter();
	int sortByIntensityAllNeighbors();
	int readStatFile();
	int readPointsList();
	int printConnectionsASCII();
	int printConnectionsASCII_trc();
	void setColorToTrackedStructures();
	void setIndexToTrackedStructures();
	void drawSegment(USrtArray& image, long time, Node* node1, Node* node2);
	virtual TimeNodes* allocTimeNodes(long time, long nbNodes) { return NULL; };
	void allocTabTimes();

public:
	int giveResult();
	TabTimes(PublicParameters* publicParameters);
	virtual ~TabTimes();
	TimeNodes* operator()(int i);
	TimeNodes** operator[](int i);
	int nbTimes;
	PublicParameters* publicParameters;

	int initialisation();
	virtual int runModel() { return T_ERR_OOP; };
	void checkLinksHomogeneity(const wchar_t* st);
};
