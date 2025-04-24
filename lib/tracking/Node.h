/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/
#pragma once

#include <cmath>

#include "ChainNeighbor.h"
#include "Parameters.h"

enum class transfer_type: uint8_t { keep_first, keep_second, give_all }; // ancien Fission/Fusion 1, 2 et 9

class Node
{
public:
	double x_, y_, z_, intensity_, mean_, surf_, color_, molId_;

	Parameters* params_;

	size_t time_;
	Node *future_, *past_;
	Node *fission1Future_, *fission2Future_, *fissionPast_;
	Node *fusion1Past_, *fusion2Past_, *fusionFuture_;
	size_t trackedIndex_;
	size_t colorIndex_[2];

	ChainNeighbor *neighborsPast_, *neighborsFuture_;
	size_t neighborsPastCount_, neighborsFutureCount_;

	//Precompute
	double distanceToBorder_, likeness_, likenessDistance_, likenessIntensity_, likenessMean_;

	Node(size_t time, Parameters* params);
	~Node();

	double Distance(const Node& node) const;
	double Likeness(const Node* n = nullptr) const;
	void UpdateConstants();
	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const;
	void ComputeBorderDistance();

	void MakeInitializationConnection();
	static void SortByIntersityNeighbors(ChainNeighbor** neighborsHead);
	double MakeChange(const double limit) { return makeChangeSimpleLinkage(limit); }

private:
	//==================== Computes ====================
	bool metropolis(double deltaE, double temperature) const;
	double sphericEffect() const;
	static double sphericEffect(const Node& n1, const Node& n2) { return n1.sphericEffect() + n2.sphericEffect(); }
	double distance(const Node& n1, const Node& n2) const;
	ChainNeighbor* getRandomNeighbor(ChainNeighbor* head, size_t size) const;

	//==================== likeness ====================
	static double getLikeness(const double difference, const double sigma) { return 1. - fabs(difference) / (3. * sigma); }
	double likenessBlink(const Node* n = nullptr) const;
	double likenessDistance(const Node* n = nullptr) const;
	double likenessIntensity(const Node* n = nullptr) const;
	double likenessMean(const Node* n = nullptr) const;

	//==================== Linkage ====================
	double makeChangeSimpleLinkage(double limit);
	bool isSimpleLinkageAvailable(const Node* node) const;
	static bool isBlinkLinkageAvailable(const Node* node);
	double simpleLinkageDisconnection(double limit) const;
	double simpleLinkageReassignment(double limit, Node* node);
	double simpleLinkageBreakAndReconnect(double limit, Node* node);
	double simpleLinkageDoubleSwap(double limit, Node* node, size_t count);
	double simpleLinkageSwap(double limit, Node* node);
};
