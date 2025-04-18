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

	virtual double Likeness(Node* n);
	virtual double LikenessBlink(Node* n);

	double DistanceToBorder() const;
	double SphericEffect() const;

	void ComputeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean, double& stdDistanceComputing,
		long& countStdDist, bool sigmaClipping, long* hist);
	double GetLikenessOfLocalNode(int lastTime);
	double MakeInitializationConnection();
	double MakeChange(double limit);


private:
	static double getLikeness(double difference, double sigma);
	double LikenessDistance(Node* n);
	double LikenessDistanceFusion(Node* n1, Node* n2);
	double LikenessFusion(Node* n1, Node* n2);
	double LikenessFission(Node* n1, Node* n2);
	double LikenessIntensity(Node* n);
	double LikenessIntensityFusion(Node* n1, Node* n2);
	double LikenessMeanIntensity(Node* n);
	double LikenessMeanIntensityFusion(Node* n1, Node* n2);
	double makeChangeSimpleLinkage(double limit);
	double makeChangeSimpleLinkage0(double limit);
	double makeChangeSimpleLinkage1(double limit, const ChainNeighbor* cRand);
	double makeChangeSimpleLinkage2(double limit, const ChainNeighbor* cRand);
	double makeChangeSimpleLinkage3(double limit, const ChainNeighbor* cRand, long count);
	double makeChangeSimpleLinkage4(double limit, const ChainNeighbor* cRand);
	double makeChangeFusion(double limit);
	double makeChangeFusion0(double limitLoss);
	double makeChangeFusion1(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFusion2(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFusion3(double limitLoss);
	double makeChangeFusion4(double limitLoss);
	double makeChangeFusion5(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFusion6(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFusion7(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFusion8(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFusion9(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission(double limitLoss);
	double makeChangeFission0(double limitLoss);
	double makeChangeFission1(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission2(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission3(double limitLoss);
	double makeChangeFission4(double limitLoss);
	double makeChangeFission5(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission6(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission7(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission8(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeFission9(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage(double limitLoss);
	double makeChangeBlinkLinkage0(double limitLoss);
	double makeChangeBlinkLinkage1(double limitLoss);
	double makeChangeBlinkLinkage2(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage3(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage4(double limitLoss, const ChainNeighbor* cRand);
	double makeChangeBlinkLinkage5(double limitLoss, const ChainNeighbor* cRand);

	static double stericEffect(const Node* n1, const Node* n2);
	static bool metropolis(double deltaE, double temperature);
};
