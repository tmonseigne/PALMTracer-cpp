/*******************************************************************************************************************
									TabTimes Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/

/********************************************************************************************************************
Class TabTimes contains several TabTimes objets for each time (nbTimes).
********************************************************************************************************************/


#include "TabTimes.h"

/********************************************************************************************************************
Constructuor(int _nbTimes) of TabTimes class. 
Create nbTimes objets TimeNodes.
********************************************************************************************************************/
TabTimes::TabTimes(Parameters* _publicParameters)
{
	parameters_ = _publicParameters;
	parameters_->nNodeAllocated_ = 0;
	parameters_->nNodeDeallocated_ = 0;
	parameters_->nChainNeighborDeallocated_ = 0;
	parameters_->nChainNeighborAllocated_ = 0;
	timeNodes_ = NULL;
	nbTimes_ = 0;
}


/********************************************************************************************************************
Destructuor() of TabTimes class. 
Delete all objets TimeNodes.
********************************************************************************************************************/
TabTimes::~TabTimes()
{
	for (int i = 0; i < nbTimes_; i++) { delete timeNodes_[i]; }
	delete [] timeNodes_;
}


/********************************************************************************************************************
allocTabTimes() alloc the timeNodes array of nbTimes long 
********************************************************************************************************************/
void TabTimes::allocTabTimes() { timeNodes_ = new TimeNodes*[nbTimes_]; }


/********************************************************************************************************************
initialisation() initialisation of the data in memory
********************************************************************************************************************/
int TabTimes::Initialisation()
{
	int ret;
	if ((ret = readPointsList()) < 0) { return ret; }

	if ((ret = linkNeighborsPresent()) < 0) { return ret; }
	if ((ret = linkNeighborsBeforeAndAfter()) < 0) { return ret; }
	if ((ret = sortByIntensityAllNeighbors()) < 0) { return ret; }

	if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1) { if ((ret = linkNeighborsBlinkBeforeAndAfter()) < 0) { return ret; } }

	CheckLinksHomogeneity(L"initalisation");
	srand(42);
	return RET_OK;
}

/********************************************************************************************************************
printConnectionsASCII_trc() print the result of the connections in a pointer as Daniel Choquet has defined:
Colonne 1         numero de la trace
Colonne 2         numero de l'image, ou temps (commencant à 1 SVP)
Colonne 3         X
Colonne 4         Y
(Colonne 5        Z)
Colonne 5 ou 6         nombre de pixels de l'objet
Colonne 6 ou 7         intensité intégrée
********************************************************************************************************************/
int TabTimes::printConnectionsASCII_trc() const
{
	unsigned long index = 0;
	for (int t = 0; t < nbTimes_; t++) { const TimeNodes* _localTimeNodes = timeNodes_[t];
		for (int n = 0; n < _localTimeNodes->nbNodes; n++) {
			Node* node = _localTimeNodes->tabNodes[n];
			if (node->connectionPast == NULL && node->trackedStructureIndex > 0) {
				for (const Node* nodeConnection = node; nodeConnection != NULL; nodeConnection = nodeConnection->connectionFuture) {
					parameters_->tracks_[index++] = node->trackedStructureIndex;
					parameters_->tracks_[index++] = nodeConnection->time + 1;
					parameters_->tracks_[index++] = nodeConnection->info->getX();
					parameters_->tracks_[index++] = nodeConnection->info->getY();

					if (Parameters::N_DIM == 3) {
						parameters_->tracks_[index++] = nodeConnection->info->getZ();
						parameters_->tracks_[index++] = nodeConnection->info->getVolume();
					}
					else { parameters_->tracks_[index++] = nodeConnection->info->getSurface(); }

					parameters_->tracks_[index++] = nodeConnection->info->getIntensity();

					parameters_->tracks_[index++] = nodeConnection->info->getZ();

					parameters_->tracks_[index++] = nodeConnection->info->getColor();
					parameters_->tracks_[index++] = nodeConnection->info->getMolId();

					// The number of trajectory
					parameters_->nTraj_ = node->trackedStructureIndex;
				}
			}
		}
	}
	return RET_OK;
}

/********************************************************************************************************************
readPointsList() reads the 'pointsList' specified in parameters_ and allocates the corresponding TimeNodes
********************************************************************************************************************/
int TabTimes::readPointsList()
{
	//unsigned long nbTimes;
	unsigned long index = 0;
	int i;

	nbTimes_ = parameters_->nTimes_;

	allocTabTimes(); //Allocation of tabTimes

	parameters_->maxX_ = 0;
	parameters_->maxY_ = 0;
	parameters_->maxZ_ = 0;

	for (int time_ = 0; time_ < (long)nbTimes_; time_++) {
		// Computing the number of nodes by frame
		unsigned long nbNodes = 0;
		for (i = index; ; /*i+=4*/ i += Parameters::N_SEG_PARAMS) {
			if ((int)parameters_->points_[i + 9] != -1 && (int)parameters_->points_[i + 9] != 0) { nbNodes++; } // surface
			else { break; }
		}


		timeNodes_[time_] = allocTimeNodes(time_, (int)nbNodes); //Allocation of timeNodes
		TimeNodes* localTimeNodes = timeNodes_[time_];

		if (localTimeNodes == NULL) { return T_ERR_ALLOC; }


		for (long node = 0; node < localTimeNodes->nbNodes; node++) {
			// Just for the case of '2D'
			InfoNode* str = (*localTimeNodes)(node)->info;
			// version correcte
			str->setX(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 3]);
			str->setY(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 4]);
			str->setZ(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 10]);

			str->setColor(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 11]);
			str->setMolId(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 12]);

			str->setIndex(node + 1);

			str->setIntensity(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 8]);

			str->setSurface((unsigned long)parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 9]);
			str->setMean(parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 8] / (unsigned long)parameters_->points_[index + node *
				             Parameters::N_SEG_PARAMS + 9]);

			if (parameters_->maxX_ < (*localTimeNodes)(node)->info->getX()) { parameters_->maxX_ = (*localTimeNodes)(node)->info->getX(); }
			if (parameters_->maxY_ < (*localTimeNodes)(node)->info->getY()) { parameters_->maxY_ = (*localTimeNodes)(node)->info->getY(); }
			if (Parameters::N_DIM == 3) { if (parameters_->maxZ_ < (*localTimeNodes)(node)->info->getZ()) { parameters_->maxZ_ = (*localTimeNodes)(node)->info->getZ(); } }
		}

		//Sorting the nodes by the intensity for one time
		int indexMax;
		for (i = 0; i < timeNodes_[time_]->nbNodes - 1; i++) {
			double maxIntensity = -1.;
			for (int j = i; j < localTimeNodes->nbNodes; j++) {
				if (maxIntensity < (*localTimeNodes)(j)->info->getIntensity()) {
					maxIntensity = (*localTimeNodes)(j)->info->getIntensity();
					indexMax = j;
				}
			}
			Node* nodeMax = localTimeNodes->tabNodes[indexMax];
			localTimeNodes->tabNodes[indexMax] = localTimeNodes->tabNodes[i];
			localTimeNodes->tabNodes[i] = nodeMax;
		}

		index += (nbNodes + 1) * Parameters::N_SEG_PARAMS; //(A vérifier)
	}
	return RET_OK;
}

/********************************************************************************************************************
operator()(int i) to access to all TimeNodes objets. 
********************************************************************************************************************/
TimeNodes* TabTimes::operator()(const int i) const { return timeNodes_[i]; }

/********************************************************************************************************************
operator[](int i) to access to pointors to all TimeNodes objets. 
********************************************************************************************************************/
TimeNodes** TabTimes::operator[](const int i) const { return timeNodes_ + i; }

/********************************************************************************************************************
int linkNeighborsPresent() to link of the neighbors in the presents of all nodes of all times
********************************************************************************************************************/
int TabTimes::linkNeighborsPresent() const
{
	int ret;
	for (long time_ = 0; time_ < nbTimes_; time_++) { if ((ret = timeNodes_[time_]->linkNeighborsPresent()) < 0) { return ret; } }
	return RET_OK;
}


/********************************************************************************************************************
int linkNeighborsBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of all times
********************************************************************************************************************/
int TabTimes::linkNeighborsBeforeAndAfter() const
{
	int ret;
	for (long time_ = 0; time_ < nbTimes_ - 1; time_++) { if ((ret = timeNodes_[time_]->linkNeighborsBeforeAndAfter(*timeNodes_[time_ + 1])) < 0) { return ret; } }
	return RET_OK;
}

/********************************************************************************************************************
int linkNeighborsBlinkBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of 
all times for blinking purpose
********************************************************************************************************************/
int TabTimes::linkNeighborsBlinkBeforeAndAfter() const
{
	int ret;
	for (long time_ = 0; time_ < nbTimes_; time_++) { for (long dt = 2; dt <= Parameters::N_TEMPORAL_WINDOWS; dt++) { if (time_ + dt < nbTimes_) { if ((ret = timeNodes_[time_]->linkNeighborsBlinkBeforeAndAfter(*timeNodes_[time_ + dt])) < 0) { return ret; } } } }
	return RET_OK;
}

/********************************************************************************************************************
int sortByIntensityAllNeighbors() sort by intensity of the node neighbors of each time 
********************************************************************************************************************/
int TabTimes::sortByIntensityAllNeighbors() const
{
	int ret;
	for (long time_ = 0; time_ < nbTimes_; time_++) { if ((ret = timeNodes_[time_]->sortByIntensityAllNeighbors()) < 0) { return ret; } }
	return RET_OK;
}


/********************************************************************************************************************
void setColorToTrackedStructures() set a color index to each tracked structures
********************************************************************************************************************/
void TabTimes::setColorToTrackedStructures() const
{
	Node* nodePastFission;
	Node* nodePastFission2;
	int newIndex = 1;
	for (long time_ = 0; time_ < nbTimes_; time_++) {
		for (long indexNode = 0; indexNode < timeNodes_[time_]->nbNodes; indexNode++) {
			Node* node = timeNodes_[time_]->tabNodes[indexNode];
			if (node->trackedStructureIndex != 0) {
				if (node->connectionPast == NULL) {
					if (node->connectionBlinkPast == NULL) {
						if (node->connectionFusion1Past != NULL) {
							if (node->likeness(node->connectionFusion1Past) > node->likeness(node->connectionFusion2Past)) {
								node->colorIndex[0] = node->connectionFusion1Past->colorIndex[0];
								node->colorIndex[1] = node->connectionFusion2Past->colorIndex[0];
							}
							else {
								node->colorIndex[0] = node->connectionFusion2Past->colorIndex[0];
								node->colorIndex[1] = node->connectionFusion1Past->colorIndex[0];
							}
						}
						else { //node->connectionFusion1Past==NULL
							if (node->connectionFissionPast != NULL) {
								Node* node2 = node == node->connectionFissionPast->connectionFission1Future ? node->connectionFissionPast->connectionFission2Future : node->connectionFissionPast->connectionFission1Future;
								for (nodePastFission = node->connectionFissionPast; nodePastFission != NULL; nodePastFission = nodePastFission->connectionPast) { nodePastFission2 = nodePastFission; }
								if (nodePastFission2->connectionFusion1Past != NULL) {
									Node* nodePastFission3 = nodePastFission2->connectionFusion1Past;
									Node* nodePastFission4 = nodePastFission2->connectionFusion2Past;
									if (node->likeness(nodePastFission3) + node2->likeness(nodePastFission4) > node2->likeness(nodePastFission3) + node->likeness(nodePastFission4)) {
										node->colorIndex[0] = nodePastFission3->colorIndex[0];
										node->colorIndex[1] = nodePastFission3->colorIndex[0];
										node2->colorIndex[0] = nodePastFission4->colorIndex[0];
										node2->colorIndex[1] = nodePastFission4->colorIndex[0];
									}
									else {
										node2->colorIndex[0] = nodePastFission3->colorIndex[0];
										node2->colorIndex[1] = nodePastFission3->colorIndex[0];
										node->colorIndex[0] = nodePastFission4->colorIndex[0];
										node->colorIndex[1] = nodePastFission4->colorIndex[0];
									}
								}
								else {
									if (node->connectionFissionPast->likeness(node) > node->connectionFissionPast->likeness(node2)) {
										node->colorIndex[0] = node->connectionFissionPast->colorIndex[0];
										node->colorIndex[1] = node->connectionFissionPast->colorIndex[0];
									}
									else {
										node->colorIndex[0] = node->connectionFissionPast->colorIndex[1];
										node->colorIndex[1] = node->connectionFissionPast->colorIndex[1];
									}
								}
							}
							else { //node->connectionFissionPast==NULL  and node->connectionFusion1Past==NULL
								if (node->connectionFission1Future != NULL) {
									node->colorIndex[0] = newIndex;
									newIndex++;
									node->colorIndex[1] = newIndex;
									newIndex++;
								}
								else {
									node->colorIndex[0] = newIndex;
									node->colorIndex[1] = newIndex;
									newIndex++;
								}
							}
						}
					}
					else {
						node->colorIndex[0] = node->connectionBlinkPast->colorIndex[0];
						node->colorIndex[1] = node->connectionBlinkPast->colorIndex[0];
					}
				}
				else { //node->connectionPast!=NULL
					node->colorIndex[0] = node->connectionPast->colorIndex[0];
					node->colorIndex[1] = node->connectionPast->colorIndex[0];
				}
				if (node->connectionFission1Future != NULL) {
					for (nodePastFission = node; nodePastFission != NULL; nodePastFission = nodePastFission->connectionPast) { nodePastFission2 = nodePastFission; }
					if (nodePastFission2->connectionFusion1Past != NULL) {
						node->colorIndex[0] = nodePastFission2->connectionFusion1Past->colorIndex[0];
						node->colorIndex[1] = nodePastFission2->connectionFusion2Past->colorIndex[0];
					}
					else {
						node->colorIndex[1] = newIndex;
						newIndex++;
					}
				}
			}
			else { //node->trackedStructureIndex==0
				node->colorIndex[0] = 0;
				node->colorIndex[1] = 0;
			}
		}
	}
}

void TabTimes::setIndexToTrackedStructures() const
{
	unsigned long newIndex = 1;
	long nbSuccessiveNodes;
	Node* node;
	Node* oldNode;
	for (long time_ = 0; time_ < nbTimes_; time_++) {
		for (long indexNode = 0; indexNode < timeNodes_[time_]->nbNodes; indexNode++) {
			if (timeNodes_[time_]->tabNodes[indexNode]->connectionPast == NULL) {
				for (node = timeNodes_[time_]->tabNodes[indexNode], oldNode = NULL, nbSuccessiveNodes = 0; node != NULL && nbSuccessiveNodes < parameters_->minLife_; node = node->connectionFuture, nbSuccessiveNodes++) {
					oldNode = node;
					node->trackedStructureIndex = 0;
				}
				bool isNewIndex = false;
				if (nbSuccessiveNodes >= parameters_->minLife_) //it was == 
				{
					isNewIndex = true;
				}
				if (timeNodes_[time_]->tabNodes[indexNode]->connectionFissionPast != NULL) { isNewIndex = true; }
				if (timeNodes_[time_]->tabNodes[indexNode]->connectionFusion1Past != NULL) { isNewIndex = true; }
				if (timeNodes_[time_]->tabNodes[indexNode]->connectionBlinkPast != NULL) { isNewIndex = true; }
				if (oldNode != NULL) {
					if (oldNode->connectionFission1Future != NULL) { isNewIndex = true; }
					if (oldNode->connectionFusionFuture != NULL) { isNewIndex = true; }
					if (oldNode->connectionBlinkFuture != NULL) { isNewIndex = true; }
				}

				if (isNewIndex) {
					timeNodes_[time_]->tabNodes[indexNode]->trackedStructureIndex = newIndex;
					newIndex++;
				}
			}
			else { timeNodes_[time_]->tabNodes[indexNode]->trackedStructureIndex = timeNodes_[time_]->tabNodes[indexNode]->connectionPast->trackedStructureIndex; }
		}
	}
}

/********************************************************************************************************************
int giveResult() just call the function to display the result
********************************************************************************************************************/
int TabTimes::GiveResult()
{
	int ret;

	setIndexToTrackedStructures();
	setColorToTrackedStructures();

	if ((ret = printConnectionsASCII_trc()) < 0) { return ret; }

	return RET_OK;
}

/********************************************************************************************************************
void checkLinksHomogeneity() just call the function to display the result
********************************************************************************************************************/
void TabTimes::CheckLinksHomogeneity(const wchar_t* st) const { for (long time_ = 0; time_ < nbTimes_; time_++) { timeNodes_[time_]->checkLinksHomogeneity(st); } }
