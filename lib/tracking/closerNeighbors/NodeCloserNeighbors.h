/*******************************************************************************************************************
									Nodes closer neighbors Class
									It is derived from the class Nodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/
#pragma once

#include "tracking/Node.h"


class NodeCloserNeighbors : public Node
{
public:
	NodeCloserNeighbors(int time, Parameters* params);

	void MakeConnection();
};
