/*******************************************************************************************************************
									TimeNodesCloserNeighbors Class
									It is derived from the class TimeNodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/

#include "TimeNodesCloserNeighbors.h"
#include "NodeCloserNeighbors.h"


/********************************************************************************************************************
allocNode() allocation of one node in the class NodeCloserNeighbors
********************************************************************************************************************/
Node* TimeNodesCloserNeighbors::allocNode() { return new NodeCloserNeighbors(time, publicParameters); }

/********************************************************************************************************************
Constructor(int _time, int _nbNodes) is run after the base constructor TimeNodes::TimeNodes
********************************************************************************************************************/
TimeNodesCloserNeighbors::TimeNodesCloserNeighbors(const int _time, const int _nbNodes, Parameters* _publicParameters): TimeNodes(_time, _nbNodes, _publicParameters) { for (int i = 0; i < nbNodes; i++) { tabNodes[i] = allocNode(); } }

/********************************************************************************************************************
makeAllConnections() finds connections from the nodes of this time to the nodes of the time after according the 
selected model
********************************************************************************************************************/
int TimeNodesCloserNeighbors::makeAllConnections() const
{
	for (int i = 0; i < nbNodes; i++) {
		NodeCloserNeighbors* node = (NodeCloserNeighbors*)tabNodes[i];
		node->makeConnection();
	}
	return RET_OK;
}
