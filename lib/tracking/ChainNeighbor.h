/*******************************************************************************************************************
									ChainNeighbor Class
								This class allows to make a chained list to chain all the presnet neighbor nodes 
*******************************************************************************************************************/
#pragma once

class Node;

class ChainNeighbor
{
public:
	ChainNeighbor(Node* neighbor, ChainNeighbor* next) : neighbor_(neighbor), next_(next) { }
	~ChainNeighbor() = default;

	Node* neighbor_;
	ChainNeighbor* next_;
};
