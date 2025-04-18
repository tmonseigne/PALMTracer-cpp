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
	ChainNeighbor(Node* _neighbor, ChainNeighbor* _next, Parameters* params);
	~ChainNeighbor();

	Parameters* parameters_;
	Node* neighbor_;
	ChainNeighbor* next_;
};
