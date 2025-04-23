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

enum class fusion_mode { from_this, from_node_to_tmp, from_node_to_this };

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

	ChainNeighbor *neighborsPast_, *neighborsFuture_;
	size_t neighborsPastCount_, neighborsFutureCount_;

	ChainNeighbor *neighborsBlinkPast_, *neighborsBlinkFuture_;
	size_t neighborsBlinkPastCount_, neighborsBlinkFutureCount_;

	Node(size_t time, Parameters* params);
	~Node();

	double Distance(const Node& node) const { return info_.Distance(node.info_); }
	double Likeness(const Node* n) const;
	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const;

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
	double makeChangeSimpleLinkage1(double limit, Node* node);
	double makeChangeSimpleLinkage2(double limit, Node* node);
	double makeChangeSimpleLinkage3(double limit, Node* node, size_t count);
	double makeChangeSimpleLinkage4(double limit, Node* node);

	double makeChangeFusion(double limit);
	double makeChangeFusionSimple(double limit, bool keepThisAsPast);
	double makeChangeFusionDouble(double limit, Node* node, bool keepFirst);
	double makeChangeFusionX(double limit, Node* node, fusion_mode mode);
	double makeChangeFusion0(double limit);
	double makeChangeFusion8(double limit, Node* node);
	double makeChangeFusion9(double limit, Node* node);

	double makeChangeFission(double limit);
	double makeChangeFissionSymmetric(double limit, Node* keepFuture);
	double makeChangeFissionReassignBoth(double limit, Node* node, Node* keep);
	double makeChangeFissionSwapFuture(double limit, Node* node, Node* swapTarget);
	double makeChangeFissionDouble(double limit, Node* node, bool keepFirst);
	double makeChangeFission0(double limit);
	double makeChangeFission9(double limit, Node* node);

	double makeChangeBlinkLinkage(double limit);
	double makeChangeBlinkLinkage0(double limit);
	double makeChangeBlinkLinkage1(double limit);
	double makeChangeBlinkLinkage2(double limit, Node* node);
	double makeChangeBlinkLinkage3(double limit, Node* node);
	double makeChangeBlinkLinkage4(double limit, Node* node);
	double makeChangeBlinkLinkage5(double limit, Node* node);

	bool metropolis(double deltaE, double temperature) const;
	double sphericEffect() const;
	static double sphericEffect(const Node& n1, const Node& n2) { return n1.sphericEffect() + n2.sphericEffect(); }
	double distance(const Node& n1, const Node& n2) const { return info_.Distance(n1.info_, n2.info_); }
	double distanceToBorder() const;
};
