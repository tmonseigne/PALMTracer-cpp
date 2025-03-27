/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/

#ifndef __nodes__
#define __nodes__

#include <cstdlib>
#include <string.h>
#include "chain_neighbor.h"
#include "general.h"
#include "info_node.h"


#define DIST_SIGNIFICANT	10.
#define DIST_INFINITY		1.e10


class Node
{
public:
	InfoNode* info;
	//int n;

	Node(int _time, PublicParameters* _publicParameters);
	~Node();

	float distance(Node* node);
	float distance(Node* n1, Node* n2);

	PublicParameters* publicParameters;
	//InfoData2D	info;

	bool checkLinksHomogeneity(const wchar_t* st);
	bool checkLinksHomogeneityInVicinity(wchar_t* st);
	bool checkLinksHomogeneityInVicinity(wchar_t* st, Node* p1, Node* p2 = NULL, Node* p3 = NULL, Node* p4 = NULL, Node* p5 = NULL, Node* p6 = NULL);
	int time;
	Node* connectionFuture;
	Node* connectionBlinkFuture;
	Node* connectionPast;
	Node* connectionBlinkPast;
	Node* connectionFission1Future;
	Node* connectionFission2Future;
	Node* connectionFissionPast;
	Node* connectionFusion1Past;
	Node* connectionFusion2Past;
	Node* connectionFusionFuture;
	unsigned long trackedStructureIndex;
	unsigned long trackedStructureAmiraIndex;
	unsigned short colorIndex[2];
	double mSD;

	ChainNeighbor* neighborsPresent;
	long neighborsPresentCount;
	ChainNeighbor* neighborsPast;
	long neighborsPastCount;
	ChainNeighbor* neighborsFuture;
	long neighborsFutureCount;

	ChainNeighbor* neighborsBlinkPast;
	long neighborsBlinkPastCount;
	ChainNeighbor* neighborsBlinkFuture;
	long neighborsBlinkFutureCount;


	int sortByIntersityNeighbors(ChainNeighbor** neighborsHead);

	virtual double likeness(Node* n) { exit(T_ERR_OOP); }
	virtual double likenessMeanIntensity(Node* n) { exit(T_ERR_OOP); }
	virtual double likenessIntensity(Node* n) { exit(T_ERR_OOP); }
	virtual double likenessDistance(Node* n) { exit(T_ERR_OOP); }
	virtual double likenessFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double likenessIntensityFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double likenessMeanIntensityFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double likenessDistanceFusion(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double likenessFission(Node* n1, Node* n2) { exit(T_ERR_OOP); }
	virtual double likenessBlink(Node* n) { exit(T_ERR_OOP); }


	double distanceToBorder();
	double stericEffect();
};


#endif//	__nodes__ 
