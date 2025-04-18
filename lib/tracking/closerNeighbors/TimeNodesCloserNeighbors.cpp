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
Node* TimeNodesCloserNeighbors::allocNode() { return new NodeCloserNeighbors(time_, parameters_); }

/********************************************************************************************************************
Constructor(int _time, int _nbNodes) is run after the base constructor TimeNodes::TimeNodes
********************************************************************************************************************/
TimeNodesCloserNeighbors::TimeNodesCloserNeighbors(const int time, const int nbNodes, Parameters* params)
	: TimeNodes(time, nbNodes, params) { for (int i = 0; i < nbNodes; i++) { nodes_[i] = allocNode(); } }

/********************************************************************************************************************
makeAllConnections() finds connections from the nodes of this time_ to the nodes of the time_ after according the selected model
********************************************************************************************************************/
int TimeNodesCloserNeighbors::MakeAllConnections() const
{
	for (int i = 0; i < nNodes_; i++) {
		NodeCloserNeighbors* node = (NodeCloserNeighbors*)nodes_[i];
		node->MakeConnection();
	}
	return RET_OK;
}
