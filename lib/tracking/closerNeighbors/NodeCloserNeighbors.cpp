/*******************************************************************************************************************
									Nodes closer neighbors Class
									It is derived from the class Nodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/


#include "NodeCloserNeighbors.h"


/********************************************************************************************************************
Constructuor(int _time) of NodeCloserNeighbors class.
********************************************************************************************************************/
NodeCloserNeighbors::NodeCloserNeighbors(const int time, Parameters* params): Node(time, params) { }

void NodeCloserNeighbors::makeConnection()
{
	Node* nMin = NULL;
	double distMin = -1.;
	for (const ChainNeighbor* c = neighborsFuture; c != NULL; c = c->next_) {
		if (c->neighbor_->connectionPast == NULL) {
			if (distMin < 0.) {
				distMin = distance(c->neighbor_);
				nMin = c->neighbor_;
			}
			else if (distance(c->neighbor_) < distMin) {
				distMin = distance(c->neighbor_);
				nMin = c->neighbor_;
			}
		}
	}
	if (nMin != NULL) {
		nMin->connectionPast = this;
		connectionFuture = nMin;
	}
}
