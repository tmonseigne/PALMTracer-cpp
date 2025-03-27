/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/

/********************************************************************************************************************
Class Node contains the caracteristic off the node, and the chain of the other nodes (time after and time before) 
links to the current objet. 
********************************************************************************************************************/

#include "nodes.h"

float globalDistMax;

/*********************************************************************************************
double distanceToBorder()
  get the distance from this to the closer border
*********************************************************************************************/
double Node::distanceToBorder()
{
	double mini = info->getX();
	mini = __min(mini, publicParameters->maxX-info->getX());
	mini = __min(mini, info->getY());
	mini = __min(mini, publicParameters->maxY-info->getY());
	return mini;
}

/*********************************************************************************************
double stericEffect()
  Compute the mean of radius of the node this
*********************************************************************************************/
double Node::stericEffect()
{
	double ret;
	double val;

	if (publicParameters->nbDim == 3) {
		val = info->getVolume();
		//volume of a shere is 4/3 * PI *r^3
		ret = pow(val * 3. / 4. / PI, 1. / 3.);
	}
	else {
		val = info->getSurface();
		//surface of a shere is PI * r^2
		ret = sqrt(val / PI);
	}
	return ret;
}

/********************************************************************************************************************
Constructuor(int _time) of Node class.
The current objets is at the time "time".
********************************************************************************************************************/
Node::Node(int _time, PublicParameters* _publicParameters)
{
	publicParameters = _publicParameters;

	time = _time;
	neighborsPresent = NULL;
	neighborsPresentCount = 0;
	neighborsPast = NULL;
	neighborsPastCount = 0;
	neighborsFuture = NULL;
	neighborsFutureCount = 0;
	neighborsBlinkPast = NULL;
	neighborsBlinkPastCount = 0;
	neighborsBlinkFuture = NULL;
	neighborsBlinkFutureCount = 0;

	connectionFuture = NULL;
	connectionPast = NULL;
	connectionBlinkPast = NULL;
	connectionBlinkFuture = NULL;
	connectionFission1Future = NULL;
	connectionFission2Future = NULL;
	connectionFissionPast = NULL;
	connectionFusion1Past = NULL;
	connectionFusion2Past = NULL;
	connectionFusionFuture = NULL;
	trackedStructureAmiraIndex = -1;
	mSD = -1;
	if (publicParameters->nbDim == 2) { info = new InfoNode2D(publicParameters); }
	else { info = new InfoNode3D(publicParameters); }

	publicParameters->nbNodeAllocated++;
}

/********************************************************************************************************************
Destructuor() of Node class. 
Delete all the linkage objets.
********************************************************************************************************************/
Node::~Node()
{
	ChainNeighbor* l;
	while (neighborsFuture != NULL) {
		l = neighborsFuture->next;
		delete neighborsFuture;
		neighborsFuture = l;
	}

	while (neighborsPast != NULL) {
		l = neighborsPast->next;
		delete neighborsPast;
		neighborsPast = l;
	}

	while (neighborsPresent != NULL) {
		l = neighborsPresent->next;
		delete neighborsPresent;
		neighborsPresent = l;
	}
	delete info;
	publicParameters->nbNodeDeallocated++;
	trackedStructureIndex = 0;
}

/********************************************************************************************************************
distance(Node *node) gives the distance between "node" and the current objet. 
********************************************************************************************************************/
float Node::distance(Node* node)
{ //Distance 1
	return info->distance(node->info);
}

/********************************************************************************************************************
distance(Node *n1, Node *n2) gives the distance between the mass center of ( n1 and n2), to the current objet. 
********************************************************************************************************************/
float Node::distance(Node* n1, Node* n2) { return info->distance(n1->info, n2->info); }


/********************************************************************************************************************
sortByIntersityNeighbors(ChainNeighbor* neighborsHead) sort the chained list begining at neighborsHead by the 
intensity. 
********************************************************************************************************************/
int Node::sortByIntersityNeighbors(ChainNeighbor** neighborsHead)
{
	bool change;
	if (*neighborsHead != NULL) {
		if ((*neighborsHead)->next != NULL) {
			do {
				change = false;
				ChainNeighbor* l = *neighborsHead;
				ChainNeighbor* lOld = NULL;
				while (l->next != NULL) {
					if (l->next->neighbor->info->getIntensity() > l->neighbor->info->getIntensity()) {
						ChainNeighbor* l2 = l->next;
						if (lOld == NULL) { *neighborsHead = l2; }
						else { lOld->next = l2; }
						l->next = l2->next;
						l2->next = l;
						change = true;
						lOld = l2;
					}
					else {
						l = l->next;
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
bool Node::checkLinksHomogeneityInVicinity(wchar_t* st)
{
	wchar_t st2[200];
	bool err = false;
	swprintf(st2, L"%ls this ", st);
	checkLinksHomogeneity(st2);
	if (connectionFuture != NULL) {
		swprintf(st2, L"%ls connectionFuture ", st);
		err = connectionFuture->checkLinksHomogeneity(st2);
	}
	if (!err && connectionFission1Future != NULL) {
		swprintf(st2, L"%ls connectionFission1Future ", st);
		err = connectionFission1Future->checkLinksHomogeneity(st2);
	}
	if (!err && connectionFission2Future != NULL) {
		swprintf(st2, L"%ls connectionFission2Future ", st);
		err = connectionFission2Future->checkLinksHomogeneity(st2);
	}
	if (!err && connectionFusionFuture != NULL) {
		swprintf(st2, L"%ls connectionFusionFuture ", st);
		err = connectionFusionFuture->checkLinksHomogeneity(st2);
	}
	if (!err && connectionBlinkFuture != NULL) {
		swprintf(st2, L"%ls connectionBlinkFuture ", st);
		err = connectionBlinkFuture->checkLinksHomogeneity(st2);
	}
	if (!err && connectionPast != NULL) {
		swprintf(st2, L"%ls connectionPast ", st);
		err = connectionPast->checkLinksHomogeneity(st2);
	}
	if (!err && connectionFusion1Past != NULL) {
		swprintf(st2, L"%ls connectionFusion1Past ", st);
		err = connectionFusion1Past->checkLinksHomogeneity(st2);
	}
	if (!err && connectionFusion2Past != NULL) {
		swprintf(st2, L"%ls connectionFusion2Past ", st);
		err = connectionFusion2Past->checkLinksHomogeneity(st2);
	}
	if (!err && connectionFissionPast != NULL) {
		swprintf(st2, L"%ls connectionFissionPast ", st);
		err = connectionFissionPast->checkLinksHomogeneity(st2);
	}
	if (!err && connectionBlinkPast != NULL) {
		swprintf(st2, L"%ls connectionBlinkPast ", st);
		err = connectionBlinkPast->checkLinksHomogeneity(st2);
	}
	return err;
}

bool Node::checkLinksHomogeneityInVicinity(wchar_t* st, Node* p1, Node* p2, Node* p3, Node* p4, Node* p5, Node* p6)
{
	wchar_t st2[256];
	bool err = false;
	if (p1 != NULL) {
		swprintf(st2, L"%ls 1 ", st);
		err = p1->checkLinksHomogeneityInVicinity(st2);
	}
	if (p2 != NULL) {
		swprintf(st2, L"%ls 2 ", st);
		err = err || p2->checkLinksHomogeneityInVicinity(st2);
	}
	if (p3 != NULL) {
		swprintf(st2, L"%ls 3 ", st);
		err = err || p3->checkLinksHomogeneityInVicinity(st2);
	}
	if (p4 != NULL) {
		swprintf(st2, L"%ls 4 ", st);
		err = err || p4->checkLinksHomogeneityInVicinity(st2);
	}
	if (p5 != NULL) {
		swprintf(st2, L"%ls 5 ", st);
		err = err || p5->checkLinksHomogeneityInVicinity(st2);
	}
	if (p6 != NULL) {
		swprintf(st2, L"%ls 6 ", st);
		err = err || p6->checkLinksHomogeneityInVicinity(st2);
	}
	return err;
}

/********************************************************************************************************************
checkLinksHomogeneity(const wchar_t* st) check all the connection to see the system homogeneity
********************************************************************************************************************/
bool Node::checkLinksHomogeneity(const wchar_t* st)
{
	bool err = false;
	if (connectionFuture != NULL) {
		if (connectionFuture->connectionPast != this) { err = true; }
		if (connectionFuture->connectionFusion1Past != NULL) { err = true; }
		if (connectionFuture->connectionFusion2Past != NULL) { err = true; }
		if (connectionFuture->connectionBlinkPast != NULL) { err = true; }
	}

	if (connectionFission1Future != NULL) {
		if (connectionFission1Future->connectionFissionPast != this) { err = true; }
		if (connectionFission1Future->connectionPast != NULL) { err = true; }
		if (connectionFission1Future->connectionFusion2Past != NULL) { err = true; }
		if (connectionFission1Future->connectionFusion1Past != NULL) { err = true; }
		if (connectionFission1Future->connectionBlinkPast != NULL) { err = true; }
	}

	if (connectionFission2Future != NULL) {
		if (connectionFission2Future->connectionFissionPast != this) { err = true; }
		if (connectionFission2Future->connectionPast != NULL) { err = true; }
		if (connectionFission2Future->connectionFusion2Past != NULL) { err = true; }
		if (connectionFission2Future->connectionFusion1Past != NULL) { err = true; }
		if (connectionFission2Future->connectionPast != NULL) { err = true; }
	}

	if (connectionFusionFuture != NULL) {
		if ((connectionFusionFuture->connectionFusion1Past != this && connectionFusionFuture->connectionFusion2Past != this) || connectionFusionFuture->connectionFusion1Past == NULL || connectionFusionFuture->connectionFusion2Past == NULL || connectionFusionFuture->connectionFusion1Past == connectionFusionFuture->connectionFusion2Past) { err = true; }
		if (connectionFusionFuture->connectionPast != NULL) { err = true; }
		if (connectionFusionFuture->connectionBlinkPast != NULL) { err = true; }
	}

	if (connectionPast != NULL) {
		if (connectionPast->connectionFuture != this) { err = true; }
		if (connectionPast->connectionFusionFuture != NULL) { err = true; }
		if (connectionPast->connectionFission1Future != NULL) { err = true; }
		if (connectionPast->connectionFission2Future != NULL) { err = true; }

		if (connectionPast->connectionFusionFuture != NULL) { err = true; }
	}

	if (connectionFissionPast != NULL) {
		if ((connectionFissionPast->connectionFission1Future != this && connectionFissionPast->connectionFission2Future != this) || connectionFissionPast->connectionFission1Future == NULL || connectionFissionPast->connectionFission2Future == NULL || connectionFissionPast->connectionFission1Future == connectionFissionPast->connectionFission2Future) { err = true; }
		if (connectionFissionPast->connectionFusionFuture != NULL) { err = true; }
		if (connectionFissionPast->connectionFuture != NULL) { err = true; }

		if (connectionFissionPast->connectionBlinkFuture != NULL) { err = true; }
	}


	if (connectionFusion1Past != NULL) {
		if (connectionFusion1Past->connectionFusionFuture != this) { err = true; }
		if (connectionFusion1Past->connectionFission1Future != NULL) { err = true; }
		if (connectionFusion1Past->connectionFission2Future != NULL) { err = true; }
		if (connectionFusion1Past->connectionFuture != NULL) { err = true; }

		if (connectionFusion1Past->connectionBlinkFuture != NULL) { err = true; }
	}

	if (connectionFusion2Past != NULL) {
		if (connectionFusion2Past->connectionFusionFuture != this) { err = true; }
		if (connectionFusion2Past->connectionFission1Future != NULL) { err = true; }
		if (connectionFusion2Past->connectionFission2Future != NULL) { err = true; }
		if (connectionFusion2Past->connectionFuture != NULL) { err = true; }

		if (connectionFusion2Past->connectionBlinkFuture != NULL) { err = true; }
	}

	if (connectionBlinkPast != NULL) {
		if (connectionBlinkPast->connectionBlinkFuture != this) { err = true; }
		if (connectionBlinkPast->connectionFission1Future != NULL) { err = true; }
		if (connectionBlinkPast->connectionFission2Future != NULL) { err = true; }
		if (connectionBlinkPast->connectionFuture != NULL) { err = true; }
		if (connectionBlinkPast->connectionFusionFuture != NULL) { err = true; }
	}

	if (connectionBlinkFuture != NULL) {
		if (connectionBlinkFuture->connectionBlinkPast != this) { err = true; }
		if (connectionBlinkFuture->connectionFusion1Past != NULL) { err = true; }
		if (connectionBlinkFuture->connectionFusion2Past != NULL) { err = true; }
		if (connectionBlinkFuture->connectionPast != NULL) { err = true; }
		if (connectionBlinkFuture->connectionFissionPast != NULL) { err = true; }
	}

	return err;
}
