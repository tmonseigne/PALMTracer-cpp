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
	Node* allocNode() const { return new Node(time_, parameters_); }

public:
	Parameters* parameters_;
	Node** nodes_;
	size_t nNodes_;
	TimeNodes(size_t time, size_t nNodes, Parameters* params);
	virtual ~TimeNodes();

	Node* operator()(const size_t i) const { return nodes_[i]; }
	void LinkNeighborsPresent() const;
	void LinkNeighborsBeforeAndAfter(const TimeNodes& futures) const;
	void LinkNeighborsBlinkBeforeAndAfter(const TimeNodes& futures) const;
	void SortByIntensityAllNeighbors() const;

	void ComputeConstantSigmas(double& stdIntensity, long& countStdInt, double& stdMean, long& countStdMean, double& stdDistance,
	                           long& countStdDist, bool sigmaClipping, long* hist = nullptr) const;
	void MakeAllInitializationConnections() const;
	double TryRandomChange(double limit) const;

	double GetLikenessTime(size_t lastTime) const;
};
