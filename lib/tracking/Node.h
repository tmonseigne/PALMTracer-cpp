/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/
#pragma once

#include "ChainNeighbor.h"
#include "InfoNode.h"
#include "Parameters.h"

class Node
{
public:
	InfoNode* info_;
	Parameters* parameters_;

	int time_;
	Node* future_, * blinkFuture_, * past_, * blinkPast_;
	Node* fission1Future_, * fission2Future_, * fissionPast_;
	Node* fusion1Past_, * fusion2Past_, * fusionFuture_;
	size_t trackedIndex_;
	size_t colorIndex_[2];

	ChainNeighbor* neighborsPresent_, * neighborsPast_, * neighborsFuture_;
	size_t neighborsPresentCount_, neighborsPastCount_, neighborsFutureCount_;

	ChainNeighbor* neighborsBlinkPast_, * neighborsBlinkFuture_;
	size_t neighborsBlinkPastCount_, neighborsBlinkFutureCount_;

	Node(int time, Parameters* params);
	virtual ~Node();

	double Distance(const Node* node) const;
	double Distance(const Node* n1, const Node* n2) const;

	bool CheckLinksHomogeneity() const;

	static int SortByIntersityNeighbors(ChainNeighbor** neighborsHead);

	virtual double Likeness(Node* n) = 0;
	virtual double LikenessMeanIntensity(Node* n) = 0;
	virtual double LikenessIntensity(Node* n) = 0;
	virtual double LikenessDistance(Node* n) = 0;
	virtual double LikenessFusion(Node* n1, Node* n2) = 0;
	virtual double LikenessIntensityFusion(Node* n1, Node* n2) = 0;
	virtual double LikenessMeanIntensityFusion(Node* n1, Node* n2) = 0;
	virtual double LikenessDistanceFusion(Node* n1, Node* n2) = 0;
	virtual double LikenessFission(Node* n1, Node* n2) = 0;
	virtual double LikenessBlink(Node* n) = 0;

	double DistanceToBorder() const;
	double SphericEffect() const;
};
