/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/
#pragma once

#include "Parameters.h"
#include "Node.h"


class TimeNodes
{
protected:
	int time_;
	virtual Node* allocNode() =0;

public:
	Parameters* parameters_;
	Node** nodes_;
	int nNodes_;
	TimeNodes(int time, int nNodes, Parameters* params);
	~TimeNodes();

	Node* operator()(int i) const;
	int LinkNeighborsPresent() const;
	int LinkNeighborsBeforeAndAfter(TimeNodes& futureTimeNodes) const;
	int LinkNeighborsBlinkBeforeAndAfter(TimeNodes& futureTimeNodes) const;
	int SortByIntensityAllNeighbors() const;

	void CheckLinksHomogeneity(const wchar_t* st) const;
};
