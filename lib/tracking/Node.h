/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/
#pragma once

#include <cstdlib>
#include <string.h>
#include "ChainNeighbor.h"
#include "InfoNode.h"
#include "Parameters.h"

class Node
{
public:
	InfoNode* info_;

	Node(int time, Parameters* params);
	~Node();

	float Distance(const Node* node) const;
	float Distance(const Node* n1, const Node* n2) const;

	Parameters* parameters_;
	//InfoData2D	info_;

	bool CheckLinksHomogeneity(const wchar_t* st) const;
	int time_;
	Node* future_;
	Node* blinkFuture_;
	Node* past_;
	Node* blinkPast_;
	Node* fission1Future_;
	Node* fission2Future_;
	Node* fissionPast_;
	Node* fusion1Past_;
	Node* fusion2Past_;
	Node* fusionFuture_;
	unsigned long trackedIndex_;
	unsigned short colorIndex_[2];

	ChainNeighbor* neighborsPresent_;
	long neighborsPresentCount_;
	ChainNeighbor* neighborsPast_;
	long neighborsPastCount_;
	ChainNeighbor* neighborsFuture_;
	long neighborsFutureCount_;

	ChainNeighbor* neighborsBlinkPast_;
	long neighborsBlinkPastCount_;
	ChainNeighbor* neighborsBlinkFuture_;
	long neighborsBlinkFutureCount_;


	static int SortByIntersityNeighbors(ChainNeighbor** neighborsHead);

	virtual double Likeness(Node* n) { exit(T_ERR_OOP); }
	virtual double LikenessMeanIntensity(Node* n) { exit(T_ERR_OOP); }
	virtual double LikenessIntensity(Node* n) { exit(T_ERR_OOP); }
	virtual double LikenessDistance(Node* n) { exit(T_ERR_OOP); }
	virtual double LikenessFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double LikenessIntensityFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double LikenessMeanIntensityFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double LikenessDistanceFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double LikenessFission(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double LikenessBlink(Node* n) { exit(T_ERR_OOP); }


	double DistanceToBorder() const;
	double StericEffect() const;
};
