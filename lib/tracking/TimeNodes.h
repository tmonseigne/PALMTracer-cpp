/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/
#pragma once

#include "Node.h"
#include "Parameters.h"


class TimeNodes
{
protected:
	size_t time_;

public:
	Parameters* params_;
	std::vector<Node> nodes_;
	size_t nNodes_;
	TimeNodes(size_t time, size_t nNodes, Parameters* params);
	virtual ~TimeNodes() { nodes_.clear(); }

	Node& operator[](const size_t i) { return nodes_[i]; }
	void LinkNeighborsBeforeAndAfter(TimeNodes& futures);
	void LinkNeighborsBlinkBeforeAndAfter(TimeNodes& futures);
	void SortByIntensityAllNeighbors();

	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const;
	void MakeAllInitializationConnections();
	double TryRandomChange(double limit);
};
