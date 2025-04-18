/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/
#pragma once

#include "Node.h"
#include "Parameters.h"


class TimeNodes
{
protected:
	int time_;
	virtual Node* allocNode() =0;

public:
	Parameters* parameters_;
	Node** nodes_;
	int nNodes_;
	TimeNodes(const int time, const int nNodes, Parameters* params) : time_(time), parameters_(params), nNodes_(nNodes) { nodes_ = new Node * [nNodes_]; }
	virtual ~TimeNodes();

	Node* operator()(int i) const;
	int LinkNeighborsPresent() const;
	int LinkNeighborsBeforeAndAfter(const TimeNodes& futures) const;
	int LinkNeighborsBlinkBeforeAndAfter(const TimeNodes& futures) const;
	int SortByIntensityAllNeighbors() const;

	void CheckLinksHomogeneity(const wchar_t* st) const;
};
