/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes 
*******************************************************************************************************************/


#include "chain_neighbor.h"

//extern Parameters *publicParameters;

/********************************************************************************************************************
ChainNeighbor(Node*	_neighbor, ChainNeighbor* _next) of ChainNeighbor class.
********************************************************************************************************************/
ChainNeighbor::ChainNeighbor(Node* _neighbor, ChainNeighbor* _next, Parameters* _publicParameters)
{
	publicParameters = _publicParameters;

	neighbor = _neighbor;
	next = _next;
	publicParameters->nChainNeighborAllocated_++;
}

/********************************************************************************************************************
Destructuor() of ChainLink class.
********************************************************************************************************************/
ChainNeighbor::~ChainNeighbor() { publicParameters->nChainNeighborDeallocated_++; }
