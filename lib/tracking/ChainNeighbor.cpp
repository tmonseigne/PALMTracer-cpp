/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes 
*******************************************************************************************************************/


#include "ChainNeighbor.h"

//extern Parameters *parameters_;

/********************************************************************************************************************
ChainNeighbor(Node*	_neighbor, ChainNeighbor* _next) of ChainNeighbor class.
********************************************************************************************************************/
ChainNeighbor::ChainNeighbor(Node* _neighbor, ChainNeighbor* _next, Parameters* params)
{
	parameters_ = params;

	neighbor_ = _neighbor;
	next_ = _next;
	parameters_->nChainNeighborAllocated_++;
}

/********************************************************************************************************************
Destructuor() of ChainLink class.
********************************************************************************************************************/
ChainNeighbor::~ChainNeighbor() { parameters_->nChainNeighborDeallocated_++; }
