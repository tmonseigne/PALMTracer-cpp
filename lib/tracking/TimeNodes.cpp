/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/

/********************************************************************************************************************
Class TimeNodes contains the Node objets for one time_.
********************************************************************************************************************/

#include "TimeNodes.h"
#include <cstdio>

/********************************************************************************************************************
Constructuor(int _time, int _nbNodes) of TimeNodes class.
The current objets is at the time_ "time_".
Create nNodes_ objets Node.
********************************************************************************************************************/
TimeNodes::TimeNodes(const int time, const int nNodes, Parameters* params)
{
	parameters_ = params;

	nNodes_ = nNodes;
	time_ = time;
	nodes_ = new Node*[nNodes_];
}


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
linkNeighbors() build the chained list of all the neighbours of each nodes presents at this time_. 
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
int TimeNodes::LinkNeighborsBeforeAndAfter(TimeNodes& futureTimeNodes) const
{
	for (long presentNode = 0; presentNode < nNodes_; presentNode++) {
		for (long futureNode = 0; futureNode < futureTimeNodes.nNodes_; futureNode++) {
			if (nodes_[presentNode]->Distance(futureTimeNodes.nodes_[futureNode]) < parameters_->maxSpeed_) {
				nodes_[presentNode]->neighborsFuture_ = new ChainNeighbor(futureTimeNodes.nodes_[futureNode], nodes_[presentNode]->neighborsFuture_,
				                                                          parameters_);
				nodes_[presentNode]->neighborsFutureCount_++;
				futureTimeNodes.nodes_[futureNode]->neighborsPast_ = new ChainNeighbor(nodes_[presentNode], futureTimeNodes.nodes_[futureNode]->neighborsPast_,
				                                                                       parameters_);
				futureTimeNodes.nodes_[futureNode]->neighborsPastCount_++;
			}
		}
	}

	return RET_OK;
}

/********************************************************************************************************************
linkBlinkBeforeAndAfter() for each structure at the time_ t, it builds the chained list of all the structures of the next 
time_ (t+dt) close to it. 
********************************************************************************************************************/
int TimeNodes::LinkNeighborsBlinkBeforeAndAfter(TimeNodes& futureTimeNodes) const
{
	for (long presentNode = 0; presentNode < nNodes_; presentNode++) {
		for (long futureNode = 0; futureNode < futureTimeNodes.nNodes_; futureNode++) {
			if (nodes_[presentNode]->Distance(futureTimeNodes.nodes_[futureNode]) < parameters_->maxSpeed_) {
				nodes_[presentNode]->neighborsBlinkFuture_ = new ChainNeighbor(futureTimeNodes.nodes_[futureNode], nodes_[presentNode]->neighborsBlinkFuture_,
				                                                               parameters_);
				nodes_[presentNode]->neighborsBlinkFutureCount_++;
				futureTimeNodes.nodes_[futureNode]->neighborsBlinkPast_ = new ChainNeighbor(nodes_[presentNode],
				                                                                            futureTimeNodes.nodes_[futureNode]->neighborsBlinkPast_,
				                                                                            parameters_);
				futureTimeNodes.nodes_[futureNode]->neighborsBlinkPastCount_++;
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
void TimeNodes::CheckLinksHomogeneity(const wchar_t* st) const { for (int i = 0; i < nNodes_; i++) { nodes_[i]->CheckLinksHomogeneity(st); } }
