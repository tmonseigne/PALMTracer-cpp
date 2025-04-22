/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/

/********************************************************************************************************************
Class TimeNodes contains the Node objets for one time.
********************************************************************************************************************/

#include "TimeNodes.h"
#include <cstdio>
#include <cstdlib>

#include "Constants.hpp"


/********************************************************************************************************************
Constructor(int _time, int _nbNodes) is run after the base constructor TimeNodes::TimeNodes
********************************************************************************************************************/
TimeNodes::TimeNodes(const size_t time, const size_t nNodes, Parameters* params)
	: time_(time), parameters_(params), nNodes_(nNodes) {
	nodes_ = new Node * [nNodes_];
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
operator()(int i) to access to all Node objets. 
********************************************************************************************************************/
Node* TimeNodes::operator()(const int i) const { return nodes_[i]; }

/********************************************************************************************************************
linkNeighbors() build the chained list of all the neighbours of each nodes presents at this time. 
********************************************************************************************************************/
int TimeNodes::LinkNeighborsPresent() const
{
	for (size_t i = 0; i < nNodes_; ++i) {
		for (size_t j = 0; j < nNodes_; ++j) {
			if (i != j) {
				if (nodes_[i]->Distance(nodes_[j]) < Parameters::DIST_MAX) {
					nodes_[i]->neighborsPresent_ = new ChainNeighbor(nodes_[j], nodes_[i]->neighborsPresent_);
					nodes_[i]->neighborsPresentCount_++;
				}
			}
		}
	}

	return RET_OK;
}


/********************************************************************************************************************
linkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+1) close to it. 
********************************************************************************************************************/
int TimeNodes::LinkNeighborsBeforeAndAfter(const TimeNodes& futures) const
{
	for (size_t present = 0; present < nNodes_; ++present) {
		for (size_t future = 0; future < futures.nNodes_; ++future) {
			if (nodes_[present]->Distance(futures.nodes_[future]) < parameters_->maxSpeed_) {
				nodes_[present]->neighborsFuture_ = new ChainNeighbor(futures.nodes_[future], nodes_[present]->neighborsFuture_);
				nodes_[present]->neighborsFutureCount_++;
				futures.nodes_[future]->neighborsPast_ = new ChainNeighbor(nodes_[present], futures.nodes_[future]->neighborsPast_);
				futures.nodes_[future]->neighborsPastCount_++;
			}
		}
	}

	return RET_OK;
}

/********************************************************************************************************************
linkBlinkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+dt) close to it. 
********************************************************************************************************************/
int TimeNodes::LinkNeighborsBlinkBeforeAndAfter(const TimeNodes& futures) const
{
	for (long present = 0; present < nNodes_; ++present) {
		for (long future = 0; future < futures.nNodes_; ++future) {
			if (nodes_[present]->Distance(futures.nodes_[future]) < parameters_->maxSpeed_) {
				nodes_[present]->neighborsBlinkFuture_ = new ChainNeighbor(futures.nodes_[future], nodes_[present]->neighborsBlinkFuture_);
				nodes_[present]->neighborsBlinkFutureCount_++;
				futures.nodes_[future]->neighborsBlinkPast_ = new ChainNeighbor(nodes_[present], futures.nodes_[future]->neighborsBlinkPast_);
				futures.nodes_[future]->neighborsBlinkPastCount_++;
			}
		}
	}

	return RET_OK;
}

/********************************************************************************************************************
sortByIntensityAllNeighbors() Sort by intersity all the neighbors of each node.
********************************************************************************************************************/
void TimeNodes::SortByIntensityAllNeighbors() const
{
	for (size_t i = 0; i < nNodes_; ++i) {
		Node::SortByIntersityNeighbors(&nodes_[i]->neighborsPast_);
		Node::SortByIntersityNeighbors(&nodes_[i]->neighborsPresent_);
		Node::SortByIntersityNeighbors(&nodes_[i]->neighborsFuture_);
	}
}


/********************************************************************************************************************
checkLinksHomogeneity() check all the connection to see the system homogeneity
********************************************************************************************************************/
void TimeNodes::CheckLinksHomogeneity() const { for (size_t i = 0; i < nNodes_; ++i) { nodes_[i]->CheckLinksHomogeneity(); } }


/********************************************************************************************************************
void computeConstantSigmas(double &stdIntensity, double &stdMean,
	double &stdDistance, long &countStd) compute the constant kMean and kIntensity
********************************************************************************************************************/
void TimeNodes::ComputeConstantSigmas(double& stdIntensity, long& countStdInt, double& stdMean, long& countStdMean,
	double& stdDistance, long& countStdDist, const bool sigmaClipping, long* hist) const
{
	for (size_t i = 0; i < nNodes_; ++i) {
		Node* node = nodes_[i];
		node->ComputeConstantSigmas(stdIntensity, countStdInt, stdMean, countStdMean, stdDistance, countStdDist, sigmaClipping, hist);
	}
}


/********************************************************************************************************************
makeAllInitializationConnections() finds connections from the nodes of this time_ to the nodes of the time_ after according to the selected model.
It returns the likeness computed for this time_
********************************************************************************************************************/
double TimeNodes::MakeAllInitializationConnections() const
{
	double globalLikeness = 0.;
	for (size_t i = 0; i < nNodes_; ++i) {
		const auto Node = nodes_[i];
		globalLikeness += Node->MakeInitializationConnection();
	}
	return globalLikeness;
}

/********************************************************************************************************************
tryRandomChange() take a random node at this time_ and try to apply a change
It returns the gain in likeness computed with this change
********************************************************************************************************************/
double TimeNodes::TryRandomChange(const double limit) const
{
	if (nNodes_ < 1) { return 0.; } //<=
	if (nNodes_ == 1) { return (nodes_[0])->MakeChange(limit); }
	return (nodes_[rand() % nNodes_])->MakeChange(limit);
}

/********************************************************************************************************************
double giveLikeness(long lastTime) compute the total likeness for all nodes
********************************************************************************************************************/
double TimeNodes::GetLikenessTime(const size_t lastTime) const
{
	double res = 0.;
	for (size_t i = 0; i < nNodes_; ++i) {
		const auto Node = nodes_[i];
		res += Node->GetLikenessOfLocalNode(lastTime);
	}
	return res;
}
