/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/

/********************************************************************************************************************
Class Node contains the caracteristic off the node, and the chain of the other nodes (time after and time before) 
links to the current objet. 
********************************************************************************************************************/

#include "Node.h"

#include <algorithm>
#include <cstdlib>
#include <numbers>

/*********************************************************************************************
double distanceToBorder()
  get the distance from this to the closer border
*********************************************************************************************/
double Node::DistanceToBorder() const
{
	const double X = info_->GetX(), Y = info_->GetY();
	const double MaxX = parameters_->maxX_, MaxY = parameters_->maxY_;
	return std::min({X, MaxX - X, Y, MaxY - Y});
}

/*********************************************************************************************
double stericEffect()
  Compute the mean of radius of the node this
*********************************************************************************************/
double Node::SphericEffect() const
{
	return sqrt(info_->GetSurface() / std::numbers::pi);
	//volume of a sphere is 4/3 * PI *r^3, surface is PI * r^2
	//if constexpr (Parameters::N_DIM == 2) { return sqrt(info_->GetSurface() / std::numbers::pi); }
	//else { return pow(info_->GetVolume() * 3.0 / 4.0 / std::numbers::pi, 1.0 / 3.0); }
}

/********************************************************************************************************************
Constructuor(int time) of Node class.
The current objets is at the time "time".
********************************************************************************************************************/
Node::Node(const int time, Parameters* params)
	: parameters_(params), time_(time),
	  future_(nullptr), blinkFuture_(nullptr), past_(nullptr), blinkPast_(nullptr),
	  fission1Future_(nullptr), fission2Future_(nullptr), fissionPast_(nullptr),
	  fusion1Past_(nullptr), fusion2Past_(nullptr), fusionFuture_(nullptr),
	  trackedIndex_(0), colorIndex_{},
	  neighborsPresent_(nullptr), neighborsPast_(nullptr), neighborsFuture_(nullptr),
	  neighborsPresentCount_(0), neighborsPastCount_(0), neighborsFutureCount_(0),
	  neighborsBlinkPast_(nullptr), neighborsBlinkFuture_(nullptr),
	  neighborsBlinkPastCount_(0), neighborsBlinkFutureCount_(0)
{
	info_ = new InfoNode2D(parameters_);
	parameters_->nNodeAllocated_++;
}

/********************************************************************************************************************
Destructuor() of Node class. 
Delete all the linkage objets.
********************************************************************************************************************/
Node::~Node()
{
	ChainNeighbor* l;
	while (neighborsFuture_ != nullptr) {
		l = neighborsFuture_->next_;
		delete neighborsFuture_;
		neighborsFuture_ = l;
	}

	while (neighborsPast_ != nullptr) {
		l = neighborsPast_->next_;
		delete neighborsPast_;
		neighborsPast_ = l;
	}

	while (neighborsPresent_ != nullptr) {
		l = neighborsPresent_->next_;
		delete neighborsPresent_;
		neighborsPresent_ = l;
	}
	delete info_;
	parameters_->nNodeDeallocated_++;
	trackedIndex_ = 0;
}

/********************************************************************************************************************
distance(Node *node) gives the distance between "node" and the current objet. 
********************************************************************************************************************/
double Node::Distance(const Node* node) const { return info_->Distance(node->info_); }

/********************************************************************************************************************
distance(Node *n1, Node *n2) gives the distance between the mass center of ( n1 and n2), to the current objet. 
********************************************************************************************************************/
double Node::Distance(const Node* n1, const Node* n2) const { return info_->Distance(n1->info_, n2->info_); }


/********************************************************************************************************************
sortByIntersityNeighbors(ChainNeighbor* neighborsHead) sort the chained list begining at neighborsHead by the 
intensity. 
********************************************************************************************************************/
int Node::SortByIntersityNeighbors(ChainNeighbor** neighborsHead)
{
	if (*neighborsHead != nullptr) {
		if ((*neighborsHead)->next_ != nullptr) {
			bool change;
			do {
				change = false;
				ChainNeighbor* l = *neighborsHead;
				ChainNeighbor* lOld = nullptr;
				while (l->next_ != nullptr) {
					if (l->next_->neighbor_->info_->GetIntensity() > l->neighbor_->info_->GetIntensity()) {
						ChainNeighbor* l2 = l->next_;
						if (lOld == nullptr) { *neighborsHead = l2; }
						else { lOld->next_ = l2; }
						l->next_ = l2->next_;
						l2->next_ = l;
						change = true;
						lOld = l2;
					}
					else {
						l = l->next_;
						if (lOld == nullptr) { lOld = *neighborsHead; }
					}
				}
			}
			while (change);
		}
	}
	return RET_OK;
}

/********************************************************************************************************************
checkLinksHomogeneity(const wchar_t* st) check all the connection to see the system homogeneity
********************************************************************************************************************/
bool Node::CheckLinksHomogeneity() const
{
	bool err = false;
	if (future_ != nullptr) {
		if (future_->past_ != this) { err = true; }
		if (future_->fusion1Past_ != nullptr) { err = true; }
		if (future_->fusion2Past_ != nullptr) { err = true; }
		if (future_->blinkPast_ != nullptr) { err = true; }
	}

	if (fission1Future_ != nullptr) {
		if (fission1Future_->fissionPast_ != this) { err = true; }
		if (fission1Future_->past_ != nullptr) { err = true; }
		if (fission1Future_->fusion2Past_ != nullptr) { err = true; }
		if (fission1Future_->fusion1Past_ != nullptr) { err = true; }
		if (fission1Future_->blinkPast_ != nullptr) { err = true; }
	}

	if (fission2Future_ != nullptr) {
		if (fission2Future_->fissionPast_ != this) { err = true; }
		if (fission2Future_->past_ != nullptr) { err = true; }
		if (fission2Future_->fusion2Past_ != nullptr) { err = true; }
		if (fission2Future_->fusion1Past_ != nullptr) { err = true; }
		if (fission2Future_->past_ != nullptr) { err = true; }
	}

	if (fusionFuture_ != nullptr) {
		if ((fusionFuture_->fusion1Past_ != this && fusionFuture_->fusion2Past_ != this) || fusionFuture_->fusion1Past_ == nullptr
		    || fusionFuture_->fusion2Past_ == nullptr || fusionFuture_->fusion1Past_ == fusionFuture_->fusion2Past_) { err = true; }
		if (fusionFuture_->past_ != nullptr) { err = true; }
		if (fusionFuture_->blinkPast_ != nullptr) { err = true; }
	}

	if (past_ != nullptr) {
		if (past_->future_ != this) { err = true; }
		if (past_->fusionFuture_ != nullptr) { err = true; }
		if (past_->fission1Future_ != nullptr) { err = true; }
		if (past_->fission2Future_ != nullptr) { err = true; }
		if (past_->fusionFuture_ != nullptr) { err = true; }
	}

	if (fissionPast_ != nullptr) {
		if ((fissionPast_->fission1Future_ != this && fissionPast_->fission2Future_ != this) || fissionPast_->fission1Future_ == nullptr
		    || fissionPast_->fission2Future_ == nullptr || fissionPast_->fission1Future_ == fissionPast_->fission2Future_) { err = true; }
		if (fissionPast_->fusionFuture_ != nullptr) { err = true; }
		if (fissionPast_->future_ != nullptr) { err = true; }
		if (fissionPast_->blinkFuture_ != nullptr) { err = true; }
	}


	if (fusion1Past_ != nullptr) {
		if (fusion1Past_->fusionFuture_ != this) { err = true; }
		if (fusion1Past_->fission1Future_ != nullptr) { err = true; }
		if (fusion1Past_->fission2Future_ != nullptr) { err = true; }
		if (fusion1Past_->future_ != nullptr) { err = true; }

		if (fusion1Past_->blinkFuture_ != nullptr) { err = true; }
	}

	if (fusion2Past_ != nullptr) {
		if (fusion2Past_->fusionFuture_ != this) { err = true; }
		if (fusion2Past_->fission1Future_ != nullptr) { err = true; }
		if (fusion2Past_->fission2Future_ != nullptr) { err = true; }
		if (fusion2Past_->future_ != nullptr) { err = true; }

		if (fusion2Past_->blinkFuture_ != nullptr) { err = true; }
	}

	if (blinkPast_ != nullptr) {
		if (blinkPast_->blinkFuture_ != this) { err = true; }
		if (blinkPast_->fission1Future_ != nullptr) { err = true; }
		if (blinkPast_->fission2Future_ != nullptr) { err = true; }
		if (blinkPast_->future_ != nullptr) { err = true; }
		if (blinkPast_->fusionFuture_ != nullptr) { err = true; }
	}

	if (blinkFuture_ != nullptr) {
		if (blinkFuture_->blinkPast_ != this) { err = true; }
		if (blinkFuture_->fusion1Past_ != nullptr) { err = true; }
		if (blinkFuture_->fusion2Past_ != nullptr) { err = true; }
		if (blinkFuture_->past_ != nullptr) { err = true; }
		if (blinkFuture_->fissionPast_ != nullptr) { err = true; }
	}

	return err;
}
