/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes 
*******************************************************************************************************************/
#pragma once

#include "Parameters.h"


class Node;

class ChainNeighbor
{
public:
	ChainNeighbor(Node* _neighbor, ChainNeighbor* _next, Parameters* _publicParameters);
	~ChainNeighbor();

	Parameters* publicParameters;
	Node* neighbor;
	ChainNeighbor* next;
};
