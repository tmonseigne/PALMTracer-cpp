/*******************************************************************************************************************
									TimeNodesDistanceIntensity Class
									It is derived from the class TimeNodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/

#include "time_nodes_dist_inten.h"


/********************************************************************************************************************
allocNode(void) allocation of one node in the class TimeNodesDistanceIntensity
********************************************************************************************************************/
Node* TimeNodesDistanceIntensity::allocNode(void) { return (new NodeDistanceIntensity(time, publicParameters)); }

/********************************************************************************************************************
Constructor(int _time, int _nbNodes) is run after the base constructor TimeNodes::TimeNodes
********************************************************************************************************************/
TimeNodesDistanceIntensity::TimeNodesDistanceIntensity(int _time, int _nbNodes, PublicParameters* _publicParameters): TimeNodes(_time, _nbNodes, _publicParameters) { for (int i = 0; i < nbNodes; i++) { tabNodes[i] = allocNode(); } }


/********************************************************************************************************************
void computeConstantSigmas(double &stdIntensityComputing, double &stdMeanComputing, 
	double &stdDistanceComputing, long &countStd) compute the constant kMean and kIntensity
********************************************************************************************************************/
void TimeNodesDistanceIntensity::computeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean, double& stdDistanceComputing, long& countStdDist, bool sigmaClipping, long* hist)
{
	NodeDistanceIntensity* node;
	for (int i = 0; i < nbNodes; i++) {
		node = (NodeDistanceIntensity*)tabNodes[i];
		node->computeConstantSigmas(stdIntensityComputing, countStdInt, stdMeanComputing, countStdMean, stdDistanceComputing, countStdDist, sigmaClipping, hist);
	}
}


/********************************************************************************************************************
makeAllInitializationConnections() finds connections from the nodes of this time to the nodes of the time after 
according the selected model.
It returns the likeness computed for this time
********************************************************************************************************************/
double TimeNodesDistanceIntensity::makeAllInitializationConnections()
{
	double globalLikeness = 0.;
	NodeDistanceIntensity* node;
	for (int i = 0; i < nbNodes; i++) {
		node = (NodeDistanceIntensity*)tabNodes[i];
		globalLikeness += node->makeInitializationConnection();
	}
	return (globalLikeness);
}

/********************************************************************************************************************
tryRandomChange() take a random node at this time and try to apply a change
It returns the gain in likeness computed with this change
********************************************************************************************************************/
double TimeNodesDistanceIntensity::tryRandomChange(double limit)
{
	if (nbNodes < 1) { return 0.; } //<=
	if (nbNodes == 1) { return (((NodeDistanceIntensity*)(tabNodes[0]))->makeChange(limit)); }
	return (((NodeDistanceIntensity*)(tabNodes[rand() % (nbNodes/*-1*/)]))->makeChange(limit));
}

/********************************************************************************************************************
double giveLikeness(long lastTime) compute the total likeness for all nodes
********************************************************************************************************************/
double TimeNodesDistanceIntensity::getLikenessTime(long lastTime)
{
	double globalLikeness = 0.;
	NodeDistanceIntensity* node;
	for (int i = 0; i < nbNodes; i++) {
		node = (NodeDistanceIntensity*)tabNodes[i];
		globalLikeness += node->getLikenessOfLocalNode(lastTime);
	}
	return (globalLikeness);
}
