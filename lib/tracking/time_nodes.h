/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/
#pragma once

#include "general.h"
#include "nodes.h"


class TimeNodes
{
protected:
	int time;
	virtual Node* allocNode() =0;

public:
	Node** tabNodes;
	TimeNodes(int _time, int _nbNodes, PublicParameters* _publicParameters);
	~TimeNodes();

	PublicParameters* publicParameters;
	Node* operator()(int i);
	int nbNodes;
	int linkNeighborsPresent();
	int linkNeighborsBeforeAndAfter(TimeNodes& futureTimeNodes);
	int linkNeighborsBlinkBeforeAndAfter(TimeNodes& futureTimeNodes);
	int sortByIntensityAllNeighbors();

	void checkLinksHomogeneity(const wchar_t* st);
};
