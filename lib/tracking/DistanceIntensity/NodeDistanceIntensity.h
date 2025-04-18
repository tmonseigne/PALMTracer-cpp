/*******************************************************************************************************************
									NodeDistanceIntensity  Class
									It is derived from the class Nodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/Node.h"
#include "tracking/Parameters.h"


class NodeDistanceIntensity final : public Node
{
private:
	static double getLikeness(double difference, double sigma);
	double LikenessDistance(Node* n) override;
	double LikenessDistanceFusion(Node* n1, Node* n2) override;
	double LikenessFusion(Node* n1, Node* n2) override;
	double LikenessFission(Node* n1, Node* n2) override;
	double LikenessIntensity(Node* n) override;
	double LikenessIntensityFusion(Node* n1, Node* n2) override;
	double LikenessMeanIntensity(Node* n) override;
	double LikenessMeanIntensityFusion(Node* n1, Node* n2) override;
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

public:
	void ComputeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean, double& stdDistanceComputing,
	                           long& countStdDist, bool sigmaClipping, long* hist);
	double GetLikenessOfLocalNode(int lastTime);
	double MakeChange(double limit);
	double Likeness(Node* n) override;
	double LikenessBlink(Node* n) override;
	NodeDistanceIntensity(int time, Parameters* params);
	double MakeInitializationConnection();
};
