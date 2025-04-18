/*******************************************************************************************************************
									Nodes closer neighbors Class
									It is derived from the class Nodes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/


#include "NodeCloserNeighbors.h"


/********************************************************************************************************************
Constructuor(int _time) of NodeCloserNeighbors class.
********************************************************************************************************************/
NodeCloserNeighbors::NodeCloserNeighbors(const int _time, Parameters* _publicParameters): Node(_time, _publicParameters) { }

void NodeCloserNeighbors::makeConnection()
{
	Node* nMin = NULL;
	double distMin = -1.;
	for (const ChainNeighbor* c = neighborsFuture; c != NULL; c = c->next) {
		if (c->neighbor->connectionPast == NULL) {
			if (distMin < 0.) {
				distMin = distance(c->neighbor);
				nMin = c->neighbor;
			}
			else if (distance(c->neighbor) < distMin) {
				distMin = distance(c->neighbor);
				nMin = c->neighbor;
			}
		}
	}
	if (nMin != NULL) {
		nMin->connectionPast = this;
		connectionFuture = nMin;
	}
}
