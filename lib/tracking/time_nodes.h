/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/
#pragma once

#include "Parameters.h"
#include "nodes.h"


class TimeNodes
{
protected:
	int time;
	virtual Node* allocNode() =0;

public:
	Node** tabNodes;
	TimeNodes(int _time, int _nbNodes, Parameters* _publicParameters);
	~TimeNodes();

	Parameters* publicParameters;
	Node* operator()(int i) const;
	int nbNodes;
	int linkNeighborsPresent() const;
	int linkNeighborsBeforeAndAfter(TimeNodes& futureTimeNodes) const;
	int linkNeighborsBlinkBeforeAndAfter(TimeNodes& futureTimeNodes) const;
	int sortByIntensityAllNeighbors() const;

	void checkLinksHomogeneity(const wchar_t* st) const;
};
