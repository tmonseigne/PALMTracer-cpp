/*******************************************************************************************************************
									TimeNodesCloserNeighbors Class
									It is derived from the class TimeNodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/
#pragma once

#include "tracking/time_nodes.h"


class TimeNodesCloserNeighbors : public TimeNodes
{
protected:
	Node* allocNode() override;

public:
	TimeNodesCloserNeighbors(int _time, int _nbNodes, PublicParameters* _publicParameters);
	int makeAllConnections();
};
