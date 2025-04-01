/*******************************************************************************************************************
									Nodes closer neighbors Class
									It is derived from the class Nodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/
#pragma once

#include "tracking/nodes.h"


class NodeCloserNeighbors : public Node
{
public:
	NodeCloserNeighbors(int _time, PublicParameters* _publicParameters);

	void makeConnection();
};

