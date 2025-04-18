/*******************************************************************************************************************
									TimeNodesCloserNeighbors Class
									It is derived from the class TimeNodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/
#pragma once

#include "tracking/TimeNodes.h"


class TimeNodesCloserNeighbors : public TimeNodes
{
protected:
	Node* allocNode() override;

public:
	TimeNodesCloserNeighbors(int time, int nbNodes, Parameters* params);
	int MakeAllConnections() const;
};
