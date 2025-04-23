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
	nodes_ = new Node*[nNodes_];
	for (size_t i = 0; i < nNodes_; ++i) { nodes_[i] = allocNode(); }
}

/********************************************************************************************************************
Destructuor() of TimeNodes class. 
Delete all objets Node.
********************************************************************************************************************/
TimeNodes::~TimeNodes()
{
	for (size_t i = 0; i < nNodes_; ++i) { delete nodes_[i]; }
	delete [] nodes_;
}

/********************************************************************************************************************
linkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+1) close to it. 
********************************************************************************************************************/
void TimeNodes::LinkNeighborsBeforeAndAfter(const TimeNodes& futures) const
{
	for (size_t present = 0; present < nNodes_; ++present) {
		for (size_t future = 0; future < futures.nNodes_; ++future) {
			if (nodes_[present]->Distance(futures.nodes_[future]) < params_->maxSpeed_) {
				nodes_[present]->neighborsFuture_ = new ChainNeighbor(futures.nodes_[future], nodes_[present]->neighborsFuture_);
				nodes_[present]->neighborsFutureCount_++;
				futures.nodes_[future]->neighborsPast_ = new ChainNeighbor(nodes_[present], futures.nodes_[future]->neighborsPast_);
				futures.nodes_[future]->neighborsPastCount_++;
			}
		}
	}
}

/********************************************************************************************************************
linkBlinkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+dt) close to it. 
********************************************************************************************************************/
void TimeNodes::LinkNeighborsBlinkBeforeAndAfter(const TimeNodes& futures) const
{
	for (size_t present = 0; present < nNodes_; ++present) {
		for (size_t future = 0; future < futures.nNodes_; ++future) {
			if (nodes_[present]->Distance(futures.nodes_[future]) < params_->maxSpeed_) {
				nodes_[present]->neighborsBlinkFuture_ = new ChainNeighbor(futures.nodes_[future], nodes_[present]->neighborsBlinkFuture_);
				nodes_[present]->neighborsBlinkFutureCount_++;
				futures.nodes_[future]->neighborsBlinkPast_ = new ChainNeighbor(nodes_[present], futures.nodes_[future]->neighborsBlinkPast_);
				futures.nodes_[future]->neighborsBlinkPastCount_++;
			}
		}
	}
}

/********************************************************************************************************************
sortByIntensityAllNeighbors() Sort by intersity all the neighbors of each node.
********************************************************************************************************************/
void TimeNodes::SortByIntensityAllNeighbors() const
{
	for (size_t i = 0; i < nNodes_; ++i) {
		Node::SortByIntersityNeighbors(&nodes_[i]->neighborsPast_);
		Node::SortByIntersityNeighbors(&nodes_[i]->neighborsFuture_);
	}
}


/********************************************************************************************************************
void computeConstantSigmas(double &intensity, double &mean,
	double &distance, long &countStd) compute the constant kMean and kIntensity
********************************************************************************************************************/
void TimeNodes::ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist,
                                      const bool clipping, long* hist) const
{
	for (size_t i = 0; i < nNodes_; ++i) { nodes_[i]->ComputeConstantSigmas(intensity, nInt, mean, nMean, distance, nDist, clipping, hist); }
}


/********************************************************************************************************************
makeAllInitializationConnections() finds connections from the nodes of this time_ to the nodes of the time_ after according to the selected model.
It returns the likeness computed for this time_
********************************************************************************************************************/
void TimeNodes::MakeAllInitializationConnections() const
{
	for (size_t i = 0; i < nNodes_; ++i) {
		const auto Node = nodes_[i];
		Node->MakeInitializationConnection();
	}
}

/********************************************************************************************************************
tryRandomChange() take a random node at this time_ and try to apply a change
It returns the gain in likeness computed with this change
********************************************************************************************************************/
double TimeNodes::TryRandomChange(const double limit) const
{
	if (nNodes_ < 1) { return 0.; } //<=
	if (nNodes_ == 1) { return nodes_[0]->MakeChange(limit); }
	return nodes_[params_->RandomModulo(nNodes_)]->MakeChange(limit);
}
