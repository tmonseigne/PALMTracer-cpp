/*******************************************************************************************************************
									Nodes closer neighbors Class
									It is derived from the class Nodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/

#ifndef __node_closer_neighbors__
#define __node_closer_neighbors__

#include "tracking/nodes.h"


class NodeCloserNeighbors : public Node
{
private:

public:
	NodeCloserNeighbors(int _time, PublicParameters* _publicParameters);

	void makeConnection(void);
};


#endif //__node_closer_neighbors__
