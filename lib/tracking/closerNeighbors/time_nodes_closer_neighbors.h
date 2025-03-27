/*******************************************************************************************************************
									TimeNodesCloserNeighbors Class
									It is derived from the class TimeNodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/

#ifndef __time_nodes_closer_neighbors__
#define __time_nodes_closer_neighbors__

#include "nodes_closer_neighbors.h"
#include "tracking/time_nodes.h"


class TimeNodesCloserNeighbors : public TimeNodes
{
protected:
	Node* allocNode(void);

public:
	TimeNodesCloserNeighbors(int _time, int _nbNodes, PublicParameters* _publicParameters);
	int makeAllConnections();
};

#endif //__time_nodes_closer_neighbors__
