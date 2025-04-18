/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/

/********************************************************************************************************************
Class Node contains the caracteristic off the node, and the chain of the other nodes (time_ after and time_ before) 
links to the current objet. 
********************************************************************************************************************/

#include "Node.h"

#include <cwchar>

#define PI 3.1415926535897932384626433832795

/*********************************************************************************************
double distanceToBorder()
  get the distance from this to the closer border
*********************************************************************************************/
double Node::DistanceToBorder() const
{
	double mini = info_->GetX();
	mini = __min(mini, parameters_->maxX_-info_->GetX());
	mini = __min(mini, info_->GetY());
	mini = __min(mini, parameters_->maxY_-info_->GetY());
	return mini;
}

/*********************************************************************************************
double stericEffect()
  Compute the mean of radius of the node this
*********************************************************************************************/
double Node::StericEffect() const
{
	double ret;
	double val;

	if (Parameters::N_DIM == 3) {
		val = info_->GetVolume();
		//volume of a shere is 4/3 * PI *r^3
		ret = pow(val * 3. / 4. / PI, 1. / 3.);
	}
	else {
		val = info_->GetSurface();
		//surface of a shere is PI * r^2
		ret = sqrt(val / PI);
	}
	return ret;
}

/********************************************************************************************************************
Constructuor(int time) of Node class.
The current objets is at the time_ "time_".
********************************************************************************************************************/
Node::Node(const int time, Parameters* params)
{
	parameters_ = params;

	time_ = time;
	neighborsPresent_ = NULL;
	neighborsPresentCount_ = 0;
	neighborsPast_ = NULL;
	neighborsPastCount_ = 0;
	neighborsFuture_ = NULL;
	neighborsFutureCount_ = 0;
	neighborsBlinkPast_ = NULL;
	neighborsBlinkPastCount_ = 0;
	neighborsBlinkFuture_ = NULL;
	neighborsBlinkFutureCount_ = 0;

	future_ = NULL;
	past_ = NULL;
	blinkPast_ = NULL;
	blinkFuture_ = NULL;
	fission1Future_ = NULL;
	fission2Future_ = NULL;
	fissionPast_ = NULL;
	fusion1Past_ = NULL;
	fusion2Past_ = NULL;
	fusionFuture_ = NULL;
	trackedAmiraIndex_ = -1;
	mSd_ = -1;
	if (Parameters::N_DIM == 2) { info_ = new InfoNode2D(parameters_); }
	else { info_ = new InfoNode3D(parameters_); }

	parameters_->nNodeAllocated_++;
}

/********************************************************************************************************************
Destructuor() of Node class. 
Delete all the linkage objets.
********************************************************************************************************************/
Node::~Node()
{
	ChainNeighbor* l;
	while (neighborsFuture_ != NULL) {
		l = neighborsFuture_->next_;
		delete neighborsFuture_;
		neighborsFuture_ = l;
	}

	while (neighborsPast_ != NULL) {
		l = neighborsPast_->next_;
		delete neighborsPast_;
		neighborsPast_ = l;
	}

	while (neighborsPresent_ != NULL) {
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
float Node::Distance(Node* node) const
{ //Distance 1
	return info_->Distance(node->info_);
}

/********************************************************************************************************************
distance(Node *n1, Node *n2) gives the distance between the mass center of ( n1 and n2), to the current objet. 
********************************************************************************************************************/
float Node::Distance(Node* n1, Node* n2) const { return info_->Distance(n1->info_, n2->info_); }


/********************************************************************************************************************
sortByIntersityNeighbors(ChainNeighbor* neighborsHead) sort the chained list begining at neighborsHead by the 
intensity. 
********************************************************************************************************************/
int Node::SortByIntersityNeighbors(ChainNeighbor** neighborsHead)
{
	bool change;
	if (*neighborsHead != NULL) {
		if ((*neighborsHead)->next_ != NULL) {
			do {
				change = false;
				ChainNeighbor* l = *neighborsHead;
				ChainNeighbor* lOld = NULL;
				while (l->next_ != NULL) {
					if (l->next_->neighbor_->info_->GetIntensity() > l->neighbor_->info_->GetIntensity()) {
						ChainNeighbor* l2 = l->next_;
						if (lOld == NULL) { *neighborsHead = l2; }
						else { lOld->next_ = l2; }
						l->next_ = l2->next_;
						l2->next_ = l;
						change = true;
						lOld = l2;
					}
					else {
						l = l->next_;
						if (lOld == NULL) { lOld = *neighborsHead; }
					}
				}
			}
			while (change);
		}
	}
	return RET_OK;
}

/********************************************************************************************************************
checkLinksHomogeneityInVicinity(wchar_t* st) check all the connection to verify the system homogeneity in the vicinity
return true if error
********************************************************************************************************************/
bool Node::CheckLinksHomogeneityInVicinity(wchar_t* st)
{
	wchar_t st2[200];
	bool err = false;
	swprintf(st2, L"%ls this ", st);
	CheckLinksHomogeneity(st2);
	if (future_ != NULL) {
		swprintf(st2, L"%ls future_ ", st);
		err = future_->CheckLinksHomogeneity(st2);
	}
	if (!err && fission1Future_ != NULL) {
		swprintf(st2, L"%ls fission1Future_ ", st);
		err = fission1Future_->CheckLinksHomogeneity(st2);
	}
	if (!err && fission2Future_ != NULL) {
		swprintf(st2, L"%ls fission2Future_ ", st);
		err = fission2Future_->CheckLinksHomogeneity(st2);
	}
	if (!err && fusionFuture_ != NULL) {
		swprintf(st2, L"%ls fusionFuture_ ", st);
		err = fusionFuture_->CheckLinksHomogeneity(st2);
	}
	if (!err && blinkFuture_ != NULL) {
		swprintf(st2, L"%ls blinkFuture_ ", st);
		err = blinkFuture_->CheckLinksHomogeneity(st2);
	}
	if (!err && past_ != NULL) {
		swprintf(st2, L"%ls past_ ", st);
		err = past_->CheckLinksHomogeneity(st2);
	}
	if (!err && fusion1Past_ != NULL) {
		swprintf(st2, L"%ls fusion1Past_ ", st);
		err = fusion1Past_->CheckLinksHomogeneity(st2);
	}
	if (!err && fusion2Past_ != NULL) {
		swprintf(st2, L"%ls fusion2Past_ ", st);
		err = fusion2Past_->CheckLinksHomogeneity(st2);
	}
	if (!err && fissionPast_ != NULL) {
		swprintf(st2, L"%ls fissionPast_ ", st);
		err = fissionPast_->CheckLinksHomogeneity(st2);
	}
	if (!err && blinkPast_ != NULL) {
		swprintf(st2, L"%ls blinkPast_ ", st);
		err = blinkPast_->CheckLinksHomogeneity(st2);
	}
	return err;
}

bool Node::CheckLinksHomogeneityInVicinity(wchar_t* st, Node* p1, Node* p2, Node* p3, Node* p4, Node* p5, Node* p6)
{
	wchar_t st2[256];
	bool err = false;
	if (p1 != NULL) {
		swprintf(st2, L"%ls 1 ", st);
		err = p1->CheckLinksHomogeneityInVicinity(st2);
	}
	if (p2 != NULL) {
		swprintf(st2, L"%ls 2 ", st);
		err = err || p2->CheckLinksHomogeneityInVicinity(st2);
	}
	if (p3 != NULL) {
		swprintf(st2, L"%ls 3 ", st);
		err = err || p3->CheckLinksHomogeneityInVicinity(st2);
	}
	if (p4 != NULL) {
		swprintf(st2, L"%ls 4 ", st);
		err = err || p4->CheckLinksHomogeneityInVicinity(st2);
	}
	if (p5 != NULL) {
		swprintf(st2, L"%ls 5 ", st);
		err = err || p5->CheckLinksHomogeneityInVicinity(st2);
	}
	if (p6 != NULL) {
		swprintf(st2, L"%ls 6 ", st);
		err = err || p6->CheckLinksHomogeneityInVicinity(st2);
	}
	return err;
}

/********************************************************************************************************************
checkLinksHomogeneity(const wchar_t* st) check all the connection to see the system homogeneity
********************************************************************************************************************/
bool Node::CheckLinksHomogeneity(const wchar_t* st) const
{
	bool err = false;
	if (future_ != NULL) {
		if (future_->past_ != this) { err = true; }
		if (future_->fusion1Past_ != NULL) { err = true; }
		if (future_->fusion2Past_ != NULL) { err = true; }
		if (future_->blinkPast_ != NULL) { err = true; }
	}

	if (fission1Future_ != NULL) {
		if (fission1Future_->fissionPast_ != this) { err = true; }
		if (fission1Future_->past_ != NULL) { err = true; }
		if (fission1Future_->fusion2Past_ != NULL) { err = true; }
		if (fission1Future_->fusion1Past_ != NULL) { err = true; }
		if (fission1Future_->blinkPast_ != NULL) { err = true; }
	}

	if (fission2Future_ != NULL) {
		if (fission2Future_->fissionPast_ != this) { err = true; }
		if (fission2Future_->past_ != NULL) { err = true; }
		if (fission2Future_->fusion2Past_ != NULL) { err = true; }
		if (fission2Future_->fusion1Past_ != NULL) { err = true; }
		if (fission2Future_->past_ != NULL) { err = true; }
	}

	if (fusionFuture_ != NULL) {
		if ((fusionFuture_->fusion1Past_ != this && fusionFuture_->fusion2Past_ != this) || fusionFuture_->fusion1Past_ == NULL || fusionFuture_->fusion2Past_
		    == NULL || fusionFuture_->fusion1Past_ == fusionFuture_->fusion2Past_) { err = true; }
		if (fusionFuture_->past_ != NULL) { err = true; }
		if (fusionFuture_->blinkPast_ != NULL) { err = true; }
	}

	if (past_ != NULL) {
		if (past_->future_ != this) { err = true; }
		if (past_->fusionFuture_ != NULL) { err = true; }
		if (past_->fission1Future_ != NULL) { err = true; }
		if (past_->fission2Future_ != NULL) { err = true; }

		if (past_->fusionFuture_ != NULL) { err = true; }
	}

	if (fissionPast_ != NULL) {
		if ((fissionPast_->fission1Future_ != this && fissionPast_->fission2Future_ != this) || fissionPast_->fission1Future_ == NULL || fissionPast_->
		    fission2Future_ == NULL || fissionPast_->fission1Future_ == fissionPast_->fission2Future_) { err = true; }
		if (fissionPast_->fusionFuture_ != NULL) { err = true; }
		if (fissionPast_->future_ != NULL) { err = true; }

		if (fissionPast_->blinkFuture_ != NULL) { err = true; }
	}


	if (fusion1Past_ != NULL) {
		if (fusion1Past_->fusionFuture_ != this) { err = true; }
		if (fusion1Past_->fission1Future_ != NULL) { err = true; }
		if (fusion1Past_->fission2Future_ != NULL) { err = true; }
		if (fusion1Past_->future_ != NULL) { err = true; }

		if (fusion1Past_->blinkFuture_ != NULL) { err = true; }
	}

	if (fusion2Past_ != NULL) {
		if (fusion2Past_->fusionFuture_ != this) { err = true; }
		if (fusion2Past_->fission1Future_ != NULL) { err = true; }
		if (fusion2Past_->fission2Future_ != NULL) { err = true; }
		if (fusion2Past_->future_ != NULL) { err = true; }

		if (fusion2Past_->blinkFuture_ != NULL) { err = true; }
	}

	if (blinkPast_ != NULL) {
		if (blinkPast_->blinkFuture_ != this) { err = true; }
		if (blinkPast_->fission1Future_ != NULL) { err = true; }
		if (blinkPast_->fission2Future_ != NULL) { err = true; }
		if (blinkPast_->future_ != NULL) { err = true; }
		if (blinkPast_->fusionFuture_ != NULL) { err = true; }
	}

	if (blinkFuture_ != NULL) {
		if (blinkFuture_->blinkPast_ != this) { err = true; }
		if (blinkFuture_->fusion1Past_ != NULL) { err = true; }
		if (blinkFuture_->fusion2Past_ != NULL) { err = true; }
		if (blinkFuture_->past_ != NULL) { err = true; }
		if (blinkFuture_->fissionPast_ != NULL) { err = true; }
	}

	return err;
}
