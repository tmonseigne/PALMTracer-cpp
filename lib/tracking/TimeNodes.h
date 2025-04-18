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
	int time_;
	virtual Node* allocNode();

public:
	Parameters* parameters_;
	Node** nodes_;
	int nNodes_;
	TimeNodes(int time, int nNodes, Parameters* params);
	virtual ~TimeNodes();

	Node* operator()(int i) const;
	int LinkNeighborsPresent() const;
	int LinkNeighborsBeforeAndAfter(const TimeNodes& futures) const;
	int LinkNeighborsBlinkBeforeAndAfter(const TimeNodes& futures) const;
	int SortByIntensityAllNeighbors() const;

	void CheckLinksHomogeneity() const;


	void ComputeConstantSigmas(double& stdIntensity, long& countStdInt, double& stdMean, long& countStdMean, double& stdDistance,
		long& countStdDist, bool sigmaClipping, long* hist = nullptr) const;
	double MakeAllInitializationConnections() const;
	double TryRandomChange(double limit) const;

	double GetLikenessTime(long lastTime) const;
};
