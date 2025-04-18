/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/

/********************************************************************************************************************
Class TimeNodes contains the Node objets for one time.
********************************************************************************************************************/

#include "TimeNodes.h"
#include <cstdio>

/********************************************************************************************************************
Destructuor() of TimeNodes class. 
Delete all objets Node.
********************************************************************************************************************/
TimeNodes::~TimeNodes()
{
	for (int i = 0; i < nNodes_; i++) { delete nodes_[i]; }
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
	for (long i = 0; i < nNodes_; i++) {
		for (long j = 0; j < nNodes_; j++) {
			if (i != j) {
				if (nodes_[i]->Distance(nodes_[j]) < Parameters::DIST_MAX) {
					nodes_[i]->neighborsPresent_ = new ChainNeighbor(nodes_[j], nodes_[i]->neighborsPresent_, parameters_);
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
	for (long present = 0; present < nNodes_; present++) {
		for (long future = 0; future < futures.nNodes_; future++) {
			if (nodes_[present]->Distance(futures.nodes_[future]) < parameters_->maxSpeed_) {
				nodes_[present]->neighborsFuture_ = new ChainNeighbor(futures.nodes_[future], nodes_[present]->neighborsFuture_, parameters_);
				nodes_[present]->neighborsFutureCount_++;
				futures.nodes_[future]->neighborsPast_ = new ChainNeighbor(nodes_[present], futures.nodes_[future]->neighborsPast_, parameters_);
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
	for (long present = 0; present < nNodes_; present++) {
		for (long future = 0; future < futures.nNodes_; future++) {
			if (nodes_[present]->Distance(futures.nodes_[future]) < parameters_->maxSpeed_) {
				nodes_[present]->neighborsBlinkFuture_ = new ChainNeighbor(futures.nodes_[future], nodes_[present]->neighborsBlinkFuture_, parameters_);
				nodes_[present]->neighborsBlinkFutureCount_++;
				futures.nodes_[future]->neighborsBlinkPast_ = new ChainNeighbor(nodes_[present], futures.nodes_[future]->neighborsBlinkPast_, parameters_);
				futures.nodes_[future]->neighborsBlinkPastCount_++;
			}
		}
	}

	return RET_OK;
}

/********************************************************************************************************************
sortByIntensityAllNeighbors() Sort by intersity all the neighbors of each nodes.
********************************************************************************************************************/
int TimeNodes::SortByIntensityAllNeighbors() const
{
	int ret;
	for (int i = 0; i < nNodes_; i++) {
		if ((ret = Node::SortByIntersityNeighbors(&nodes_[i]->neighborsPast_)) < 0) { return ret; }
		if ((ret = Node::SortByIntersityNeighbors(&nodes_[i]->neighborsPresent_)) < 0) { return ret; }
		if ((ret = Node::SortByIntersityNeighbors(&nodes_[i]->neighborsFuture_)) < 0) { return ret; }
	}
	return RET_OK;
}


/********************************************************************************************************************
checkLinksHomogeneity() check all the connection to see the system homogeneity
********************************************************************************************************************/
void TimeNodes::CheckLinksHomogeneity() const { for (int i = 0; i < nNodes_; i++) { nodes_[i]->CheckLinksHomogeneity(); } }
