/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/
#pragma once

#include <array>
#include <cmath>

#include "ChainNeighbor.h"
#include "InfoNode.h"
#include "Parameters.h"

class Node
{
public:
	InfoNode info_;
	Parameters* params_;

	size_t time_;
	Node *future_, *blinkFuture_, *past_, *blinkPast_;
	Node *fission1Future_, *fission2Future_, *fissionPast_;
	Node *fusion1Past_, *fusion2Past_, *fusionFuture_;
	size_t trackedIndex_;
	size_t colorIndex_[2];

	ChainNeighbor *neighborsPresent_, *neighborsPast_, *neighborsFuture_;
	size_t neighborsPresentCount_, neighborsPastCount_, neighborsFutureCount_;

	ChainNeighbor *neighborsBlinkPast_, *neighborsBlinkFuture_;
	size_t neighborsBlinkPastCount_, neighborsBlinkFutureCount_;

	Node(size_t time, Parameters* params);
	~Node();

	double Distance(const Node* node) const { return info_.Distance(node->info_); }
	double Likeness(const Node* n) const;
	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance,
	                           size_t& nDist, bool clipping, long* hist) const;

	void MakeInitializationConnection();
	static void SortByIntersityNeighbors(ChainNeighbor** neighborsHead);
	double MakeChange(double limit);

private:
	static double getLikeness(const double difference, const double sigma) { return 1. - fabs(difference) / (3. * sigma); }
	double likenessBlink(const Node* n) const;
	double likenessDistance(const Node* n) const;
	double likenessDistanceFusion(const Node* n1, const Node* n2) const;
	double likenessFusion(const Node* n1, const Node* n2) const;
	double likenessIntensity(const Node* n) const;
	double likenessIntensityFusion(const Node* n1, const Node* n2) const;
	double likenessMeanIntensity(const Node* n) const;
	double likenessMeanIntensityFusion(const Node* n1, const Node* n2) const;

	double makeChangeSimpleLinkage(double limit);
	double makeChangeSimpleLinkage0(double limit) const;
	double makeChangeSimpleLinkage1(double limit, const ChainNeighbor* cRand);
	double makeChangeSimpleLinkage2(double limit, const ChainNeighbor* cRand);
	double makeChangeSimpleLinkage3(double limit, const ChainNeighbor* cRand, size_t count);
	double makeChangeSimpleLinkage4(double limit, const ChainNeighbor* cRand);

	double makeChangeFusion(double limit);
	double makeChangeFusion0(double limit);
	double makeChangeFusion1(double limit, const ChainNeighbor* cRand);
	double makeChangeFusion2(double limit, const ChainNeighbor* cRand);
	double makeChangeFusion3(double limit);
	double makeChangeFusion4(double limit);
	double makeChangeFusion5(double limit, const ChainNeighbor* cRand);
	double makeChangeFusion6(double limit, const ChainNeighbor* cRand) const;
	double makeChangeFusion7(double limit, const ChainNeighbor* cRand);
	double makeChangeFusion8(double limit, const ChainNeighbor* cRand);
	double makeChangeFusion9(double limit, const ChainNeighbor* cRand);

	double makeChangeFission(double limit);
	double makeChangeFission0(double limit);
	double makeChangeFission1(double limit, const ChainNeighbor* cRand);
	double makeChangeFission2(double limit, const ChainNeighbor* cRand);
	double makeChangeFission3(double limit);
	double makeChangeFission4(double limit);
	double makeChangeFission5(double limit, const ChainNeighbor* cRand);
	double makeChangeFission6(double limit, const ChainNeighbor* cRand);
	double makeChangeFission7(double limit, const ChainNeighbor* cRand);
	double makeChangeFission8(double limit, const ChainNeighbor* cRand);
	double makeChangeFission9(double limit, const ChainNeighbor* cRand);

	double makeChangeBlinkLinkage(double limit);
	double makeChangeBlinkLinkage0(double limit);
	double makeChangeBlinkLinkage1(double limit);
	double makeChangeBlinkLinkage2(double limit, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage3(double limit, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage4(double limit, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage5(double limit, const ChainNeighbor* cRand);

	bool metropolis(double deltaE, double temperature) const;
	double sphericEffect() const;
	static double sphericEffect(const Node* n1, const Node* n2) { return n1->sphericEffect() + n2->sphericEffect(); }
	double distance(const Node* n1, const Node* n2) const { return info_.Distance(n1->info_, n2->info_); }
	double distanceToBorder() const;
};
