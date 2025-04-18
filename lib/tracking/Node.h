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


#define DIST_SIGNIFICANT	10.
#define DIST_INFINITY		1.e10


class Node
{
public:
	InfoNode* info_;

	Node(int time, Parameters* params);
	~Node();

	float Distance(Node* node) const;
	float Distance(Node* n1, Node* n2) const;

	Parameters* parameters_;
	//InfoData2D	info_;

	bool CheckLinksHomogeneity(const wchar_t* st) const;
	bool CheckLinksHomogeneityInVicinity(wchar_t* st);
	static bool CheckLinksHomogeneityInVicinity(wchar_t* st, Node* p1, Node* p2 = nullptr, Node* p3 = nullptr, Node* p4 = nullptr, Node* p5 = nullptr, Node* p6 = nullptr);
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
	unsigned long trackedAmiraIndex_;
	unsigned short colorIndex_[2];
	double mSd_;

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
