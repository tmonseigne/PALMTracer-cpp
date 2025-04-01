/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes 
*******************************************************************************************************************/
#pragma once

#include "general.h"


class Node;

class ChainNeighbor
{
public:
	ChainNeighbor(Node* _neighbor, ChainNeighbor* _next, PublicParameters* _publicParameters);
	~ChainNeighbor();

	PublicParameters* publicParameters;
	Node* neighbor;
	ChainNeighbor* next;
};
