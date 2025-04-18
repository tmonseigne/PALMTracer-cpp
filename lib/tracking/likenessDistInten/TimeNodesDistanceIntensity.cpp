/*******************************************************************************************************************
									TimeNodesDistanceIntensity Class
									It is derived from the class TimeNodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/

#include "NodeDistanceIntensity.h"
#include "TimeNodesDistanceIntensity.h"


/********************************************************************************************************************
allocNode() allocation of one node in the class TimeNodesDistanceIntensity
********************************************************************************************************************/
Node* TimeNodesDistanceIntensity::allocNode() { return new NodeDistanceIntensity(time_, parameters_); }

/********************************************************************************************************************
Constructor(int _time, int _nbNodes) is run after the base constructor TimeNodes::TimeNodes
********************************************************************************************************************/
TimeNodesDistanceIntensity::TimeNodesDistanceIntensity(const int time, const int nNodes, Parameters* params)
	: TimeNodes(time, nNodes, params) { for (int i = 0; i < nNodes_; i++) { nodes_[i] = allocNode(); } }


/********************************************************************************************************************
void computeConstantSigmas(double &stdIntensityComputing, double &stdMeanComputing, 
	double &stdDistanceComputing, long &countStd) compute the constant kMean and kIntensity
********************************************************************************************************************/
void TimeNodesDistanceIntensity::ComputeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean,
                                                       double& stdDistanceComputing, long& countStdDist, const bool sigmaClipping, long* hist) const
{
	for (int i = 0; i < nNodes_; i++) {
		NodeDistanceIntensity* node = (NodeDistanceIntensity*)nodes_[i];
		node->ComputeConstantSigmas(stdIntensityComputing, countStdInt, stdMeanComputing, countStdMean, stdDistanceComputing, countStdDist, sigmaClipping,
		                            hist);
	}
}


/********************************************************************************************************************
makeAllInitializationConnections() finds connections from the nodes of this time_ to the nodes of the time_ after 
according the selected model.
It returns the likeness computed for this time_
********************************************************************************************************************/
double TimeNodesDistanceIntensity::MakeAllInitializationConnections() const
{
	double globalLikeness = 0.;
	for (int i = 0; i < nNodes_; i++) {
		NodeDistanceIntensity* node = (NodeDistanceIntensity*)nodes_[i];
		globalLikeness += node->MakeInitializationConnection();
	}
	return globalLikeness;
}

/********************************************************************************************************************
tryRandomChange() take a random node at this time_ and try to apply a change
It returns the gain in likeness computed with this change
********************************************************************************************************************/
double TimeNodesDistanceIntensity::TryRandomChange(const double limit) const
{
	if (nNodes_ < 1) { return 0.; } //<=
	if (nNodes_ == 1) { return ((NodeDistanceIntensity*)nodes_[0])->MakeChange(limit); }
	return ((NodeDistanceIntensity*)nodes_[rand() % nNodes_/*-1*/])->MakeChange(limit);
}

/********************************************************************************************************************
double giveLikeness(long lastTime) compute the total likeness for all nodes
********************************************************************************************************************/
double TimeNodesDistanceIntensity::GetLikenessTime(const long lastTime) const
{
	double globalLikeness = 0.;
	for (int i = 0; i < nNodes_; i++) {
		NodeDistanceIntensity* node = (NodeDistanceIntensity*)nodes_[i];
		globalLikeness += node->GetLikenessOfLocalNode(lastTime);
	}
	return globalLikeness;
}
