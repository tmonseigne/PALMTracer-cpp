#include "stdafx.h"


/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes 
*******************************************************************************************************************/


#include "chain_neighbor.h"

//extern PublicParameters *publicParameters;

/********************************************************************************************************************
ChainNeighbor(Node*	_neighbor, ChainNeighbor* _next) of ChainNeighbor class.
********************************************************************************************************************/
ChainNeighbor::ChainNeighbor(Node* _neighbor, ChainNeighbor* _next, PublicParameters* _publicParameters)
{
	publicParameters = _publicParameters;

	neighbor = _neighbor;
	next = _next;
	publicParameters->nbChainNeighborAllocated++;
}

/********************************************************************************************************************
Destructuor() of ChainLink class.
********************************************************************************************************************/
ChainNeighbor::~ChainNeighbor() { publicParameters->nbChainNeighborDeallocated++; }
