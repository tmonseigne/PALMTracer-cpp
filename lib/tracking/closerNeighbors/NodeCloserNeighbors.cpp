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

void NodeCloserNeighbors::MakeConnection()
{
	Node* nMin = nullptr;
	double distMin = -1.;
	for (const ChainNeighbor* c = neighborsFuture_; c != nullptr; c = c->next_) {
		if (c->neighbor_->past_ == nullptr) {
			if (distMin < 0.) {
				distMin = Distance(c->neighbor_);
				nMin = c->neighbor_;
			}
			else if (Distance(c->neighbor_) < distMin) {
				distMin = Distance(c->neighbor_);
				nMin = c->neighbor_;
			}
		}
	}
	if (nMin != nullptr) {
		nMin->past_ = this;
		future_ = nMin;
	}
}
