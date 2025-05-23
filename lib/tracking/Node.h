/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/
#pragma once

#include <array>
#include <cmath>

#include "Parameters.hpp"

namespace palm_tracer {
class Node {
public:
	double x_, y_, z_, intensity_, mean_, surf_, color_, molId_;

	Parameters* params_;

	size_t time_;
	Node *future_, *past_;
	size_t trackedIndex_;
	size_t colorIndex_[2];

	std::array<std::vector<Node*>, 2> neighbors_;

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
	double MakeChange(const double limit) { return makeChangeSimpleLinkage(limit); }

private:
	//==================== Computes ====================
	bool metropolis(double deltaE, double temperature) const;
	double sphericEffect() const;
	static double sphericEffect(const Node& n1, const Node& n2) { return n1.sphericEffect() + n2.sphericEffect(); }
	double distance(const Node& n1, const Node& n2) const;
	Node* getRandomNeighbor(const size_t type = 0) const { return neighbors_[type][params_->RandomModulo(neighbors_[type].size())]; }

	//==================== likeness ====================
	static double getLikeness(const double difference, const double sigma) { return 1. - std::fabs(difference) / (3. * sigma); }
	double likenessBlink(const Node* n = nullptr) const;
	double likenessDistance(const Node* n = nullptr) const;
	double likenessIntensity(const Node* n = nullptr) const;
	double likenessMean(const Node* n = nullptr) const;

	//==================== Linkage ====================
	double makeChangeSimpleLinkage(double limit);
	bool isSimpleLinkageAvailable(const Node* node) const;
	double simpleLinkageDisconnection(double limit) const;
	double simpleLinkageReassignment(double limit, Node* node);
	double simpleLinkageBreakAndReconnect(double limit, Node* node);
	double simpleLinkageSwap(double limit, Node* node);
	double simpleLinkageDoubleSwap(double limit, Node* node, size_t count);
};
}
