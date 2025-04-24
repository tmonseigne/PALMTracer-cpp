/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/

/********************************************************************************************************************
Class TimeNodes contains the Node objets for one time.
********************************************************************************************************************/

#include "TimeNodes.h"

/********************************************************************************************************************
Constructor(int _time, int _nbNodes) is run after the base constructor TimeNodes::TimeNodes
********************************************************************************************************************/
TimeNodes::TimeNodes(const size_t time, const size_t nNodes, Parameters* params)
	: time_(time), params_(params), nNodes_(nNodes)
{
	nodes_.reserve(nNodes_);
	for (size_t i = 0; i < nNodes_; ++i) { nodes_.emplace_back(time_, params_); }
}

/********************************************************************************************************************
linkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+1) close to it. 
********************************************************************************************************************/
void TimeNodes::LinkNeighborsBeforeAndAfter(TimeNodes& futures)
{
	for (auto& present : nodes_) {
		for (auto& future : futures.nodes_) {
			if (present.Distance(future) < params_->maxSpeed_) {
				present.neighborsFuture_ = new ChainNeighbor(&future, present.neighborsFuture_);
				present.neighborsFutureCount_++;
				future.neighborsPast_ = new ChainNeighbor(&present, future.neighborsPast_);
				future.neighborsPastCount_++;
			}
		}
	}
}

/********************************************************************************************************************
linkBlinkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+dt) close to it. 
********************************************************************************************************************/
void TimeNodes::LinkNeighborsBlinkBeforeAndAfter(TimeNodes& futures)
{
	for (auto& present : nodes_) {
		for (auto& future : futures.nodes_) {
			if (present.Distance(future) < params_->maxSpeed_) {
				present.neighborsBlinkFuture_ = new ChainNeighbor(&future, present.neighborsBlinkFuture_);
				present.neighborsBlinkFutureCount_++;
				future.neighborsBlinkPast_ = new ChainNeighbor(&present, future.neighborsBlinkPast_);
				future.neighborsBlinkPastCount_++;
			}
		}
	}
}

/********************************************************************************************************************
sortByIntensityAllNeighbors() Sort by intersity all the neighbors of each node.
********************************************************************************************************************/
void TimeNodes::SortByIntensityAllNeighbors()
{
	for (auto& n : nodes_) {
		Node::SortByIntersityNeighbors(&n.neighborsPast_);
		Node::SortByIntersityNeighbors(&n.neighborsFuture_);
	}
}


/********************************************************************************************************************
void computeConstantSigmas(double &intensity, double &mean,
	double &distance, long &countStd) compute the constant kMean and kIntensity
********************************************************************************************************************/
void TimeNodes::ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const
{
	for (const auto& N : nodes_) { N.ComputeConstantSigmas(intensity, nInt, mean, nMean, distance, nDist); }
}
