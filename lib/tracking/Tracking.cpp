/*******************************************************************************************************************
									Tracking Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/

/********************************************************************************************************************
Class Tracking contains several Tracking objets for each time (nbTimes).
********************************************************************************************************************/


#include "Tracking.h"
#include <algorithm>

#include "Constants.hpp"

/********************************************************************************************************************
Destructuor() of Tracking class. 
Delete all objets TimeNodes.
********************************************************************************************************************/
Tracking::~Tracking()
{
	for (size_t i = 0; i < nbTimes_; ++i) { delete timeNodes_[i]; }
	delete [] timeNodes_;
}

/********************************************************************************************************************
allocTabTimes() alloc the timeNodes array of nbTimes long 
********************************************************************************************************************/
void Tracking::allocTabTimes() { timeNodes_ = new TimeNodes*[nbTimes_]; }

/********************************************************************************************************************
initialisation() initialisation of the data in memory
********************************************************************************************************************/
void Tracking::Initialisation()
{
	readPointsList();
	linkNeighborsPresent();
	linkNeighborsBeforeAndAfter();
	sortByIntensityAllNeighbors();
	if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1) { linkNeighborsBlinkBeforeAndAfter(); }
	srand(42);
}

int Tracking::RunModel()
{
	{
		double error;

		constexpr long CountLearning = 5;

		for (long iLearning = 0; iLearning < CountLearning; iLearning++) {
			if (iLearning == CountLearning - 1) { error = 20; }
			else { error = 0; }
			long cpt = 0;

			if (iLearning == 0) {
				computeConstantes();
				for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t]->MakeAllInitializationConnections(); }
			}
			else { updateConstantes(); }

			double totalLikeness = getLikenessGlobal();

			int nbInt = 0;

			int bcl = 0;
			do {
				if (iLearning == CountLearning - 1) { error *= 1. - parameters_->decrease_ / 100.; }

				nbInt++;

				double deltaLikeness = 0;
				for (size_t i = 0; i < Parameters::N_ITERATIONS; ++i) {
					for (size_t t = 0; t < nbTimes_ - 1; ++t) { deltaLikeness += timeNodes_[t]->TryRandomChange(error); }
				}
				totalLikeness += deltaLikeness;

				if (fabs(deltaLikeness) < 1e-8) { cpt++; }
				else { cpt = 0; }

				bcl++;
			}
			while (cpt < 10);
		}

		return RET_OK;
	}
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
void Tracking::printConnectionsASCII_trc() const
{
	size_t index = 0;
	for (size_t t = 0; t < nbTimes_; ++t) {
		const TimeNodes* localTimeNodes = timeNodes_[t];
		for (size_t n = 0; n < localTimeNodes->nNodes_; n++) {
			const Node* node = localTimeNodes->nodes_[n];
			if (node->past_ == nullptr && node->trackedIndex_ > 0) {
				for (const Node* nodeConnection = node; nodeConnection != nullptr; nodeConnection = nodeConnection->future_) {
					parameters_->tracks_[index++] = static_cast<double>(node->trackedIndex_);
					parameters_->tracks_[index++] = static_cast<double>(nodeConnection->time_ + 1);
					parameters_->tracks_[index++] = nodeConnection->info_.x_;
					parameters_->tracks_[index++] = nodeConnection->info_.y_;
					parameters_->tracks_[index++] = nodeConnection->info_.surf_;
					parameters_->tracks_[index++] = nodeConnection->info_.intensity_;
					parameters_->tracks_[index++] = nodeConnection->info_.z_;
					parameters_->tracks_[index++] = nodeConnection->info_.color_;
					parameters_->tracks_[index++] = nodeConnection->info_.molId_;

					// The number of trajectory
					parameters_->nTraj_ = node->trackedIndex_;
				}
			}
		}
	}
}

/********************************************************************************************************************
readPointsList() reads the 'pointsList' specified in parameters_ and allocates the corresponding TimeNodes
********************************************************************************************************************/
void Tracking::readPointsList()
{
	allocTabTimes(); //Allocation of tabTimes
	parameters_->maxX_ = 0;
	parameters_->maxY_ = 0;
	parameters_->maxZ_ = 0;

	for (size_t t = 0, index = 0; t < nbTimes_; t++) {
		// Computing the number of nodes by frame
		size_t nbNodes = 0;
		for (size_t i = index; ; i += Parameters::N_SEG_PARAMS) {
			if (static_cast<int>(parameters_->points_[i + 9]) != -1 && static_cast<int>(parameters_->points_[i + 9]) != 0) { nbNodes++; } // surface
			else { break; }
		}


		timeNodes_[t] = allocTimeNodes(t, nbNodes); //Allocation of timeNodes
		const TimeNodes* localTimeNodes = timeNodes_[t];

		for (size_t node = 0; node < localTimeNodes->nNodes_; node++) {
			// Just for the case of '2D'
			InfoNode& str = (*localTimeNodes)(node)->info_;
			// version correcte
			str.x_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 3];
			str.y_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 4];
			str.z_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 10];

			str.color_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 11];
			str.molId_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 12];

			str.intensity_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 8];

			str.surf_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 9];
			str.mean_ = parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 8] / parameters_->points_[index + node * Parameters::N_SEG_PARAMS + 9];

			parameters_->maxX_ = std::max(parameters_->maxX_, (*localTimeNodes)(node)->info_.x_);
			parameters_->maxY_ = std::max(parameters_->maxY_, (*localTimeNodes)(node)->info_.y_);
		}

		//Sorting the nodes by the intensity for one time
		size_t indexMax = 0;
		for (size_t i = 0; i < timeNodes_[t]->nNodes_ - 1; ++i) {
			double maxIntensity = -1.;
			for (size_t j = i; j < localTimeNodes->nNodes_; ++j) {
				if (maxIntensity < (*localTimeNodes)(j)->info_.intensity_) {
					maxIntensity = (*localTimeNodes)(j)->info_.intensity_;
					indexMax = j;
				}
			}
			Node* nodeMax = localTimeNodes->nodes_[indexMax];
			localTimeNodes->nodes_[indexMax] = localTimeNodes->nodes_[i];
			localTimeNodes->nodes_[i] = nodeMax;
		}

		index += (nbNodes + 1) * Parameters::N_SEG_PARAMS; //(A vérifier)
	}
}

/********************************************************************************************************************
void linkNeighborsBlinkBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of 
all times for blinking purpose
********************************************************************************************************************/
void Tracking::linkNeighborsBlinkBeforeAndAfter() const
{
	for (size_t t = 0; t < nbTimes_; ++t) {
		for (size_t dt = 2; dt <= Parameters::N_TEMPORAL_WINDOWS; dt++) {
			if (t + dt < nbTimes_) { timeNodes_[t]->LinkNeighborsBlinkBeforeAndAfter(*timeNodes_[t + dt]); }
		}
	}
}

/********************************************************************************************************************
int sortByIntensityAllNeighbors() sort by intensity of the node neighbors of each time 
********************************************************************************************************************/
void Tracking::sortByIntensityAllNeighbors() const { for (size_t t = 0; t < nbTimes_; ++t) { timeNodes_[t]->SortByIntensityAllNeighbors(); } }


/********************************************************************************************************************
void setColorToTrackedStructures() set a color index to each tracked structures
********************************************************************************************************************/
void Tracking::setColorToTrackedStructures() const
{
	Node* nodePastFission;
	Node* nodePastFission2;
	int newIndex = 1;
	for (size_t time = 0; time < nbTimes_; time++) {
		for (size_t indexNode = 0; indexNode < timeNodes_[time]->nNodes_; indexNode++) {
			Node* node = timeNodes_[time]->nodes_[indexNode];
			if (node->trackedIndex_ != 0) {
				if (node->past_ == nullptr) {
					if (node->blinkPast_ == nullptr) {
						if (node->fusion1Past_ != nullptr) {
							if (node->Likeness(node->fusion1Past_) > node->Likeness(node->fusion2Past_)) {
								node->colorIndex_[0] = node->fusion1Past_->colorIndex_[0];
								node->colorIndex_[1] = node->fusion2Past_->colorIndex_[0];
							}
							else {
								node->colorIndex_[0] = node->fusion2Past_->colorIndex_[0];
								node->colorIndex_[1] = node->fusion1Past_->colorIndex_[0];
							}
						}
						else { //node->fusion1Past_==NULL
							if (node->fissionPast_ != nullptr) {
								Node* node2 = node == node->fissionPast_->fission1Future_
									              ? node->fissionPast_->fission2Future_
									              : node->fissionPast_->fission1Future_;

								for (nodePastFission = node->fissionPast_; nodePastFission != nullptr; nodePastFission = nodePastFission->past_) {
									nodePastFission2 = nodePastFission;
								}
								if (nodePastFission2->fusion1Past_ != nullptr) {
									Node* nodePastFission3 = nodePastFission2->fusion1Past_;
									Node* nodePastFission4 = nodePastFission2->fusion2Past_;
									if (node->Likeness(nodePastFission3) + node2->Likeness(nodePastFission4) > node2->Likeness(nodePastFission3) + node->
									    Likeness(nodePastFission4)) {
										node->colorIndex_[0] = nodePastFission3->colorIndex_[0];
										node->colorIndex_[1] = nodePastFission3->colorIndex_[0];
										node2->colorIndex_[0] = nodePastFission4->colorIndex_[0];
										node2->colorIndex_[1] = nodePastFission4->colorIndex_[0];
									}
									else {
										node2->colorIndex_[0] = nodePastFission3->colorIndex_[0];
										node2->colorIndex_[1] = nodePastFission3->colorIndex_[0];
										node->colorIndex_[0] = nodePastFission4->colorIndex_[0];
										node->colorIndex_[1] = nodePastFission4->colorIndex_[0];
									}
								}
								else {
									if (node->fissionPast_->Likeness(node) > node->fissionPast_->Likeness(node2)) {
										node->colorIndex_[0] = node->fissionPast_->colorIndex_[0];
										node->colorIndex_[1] = node->fissionPast_->colorIndex_[0];
									}
									else {
										node->colorIndex_[0] = node->fissionPast_->colorIndex_[1];
										node->colorIndex_[1] = node->fissionPast_->colorIndex_[1];
									}
								}
							}
							else { //node->fissionPast_==NULL  and node->fusion1Past_==NULL
								if (node->fission1Future_ != nullptr) {
									node->colorIndex_[0] = newIndex;
									newIndex++;
									node->colorIndex_[1] = newIndex;
									newIndex++;
								}
								else {
									node->colorIndex_[0] = newIndex;
									node->colorIndex_[1] = newIndex;
									newIndex++;
								}
							}
						}
					}
					else {
						node->colorIndex_[0] = node->blinkPast_->colorIndex_[0];
						node->colorIndex_[1] = node->blinkPast_->colorIndex_[0];
					}
				}
				else { //node->past_!=NULL
					node->colorIndex_[0] = node->past_->colorIndex_[0];
					node->colorIndex_[1] = node->past_->colorIndex_[0];
				}
				if (node->fission1Future_ != nullptr) {
					for (nodePastFission = node; nodePastFission != nullptr; nodePastFission = nodePastFission->past_) { nodePastFission2 = nodePastFission; }
					if (nodePastFission2->fusion1Past_ != nullptr) {
						node->colorIndex_[0] = nodePastFission2->fusion1Past_->colorIndex_[0];
						node->colorIndex_[1] = nodePastFission2->fusion2Past_->colorIndex_[0];
					}
					else {
						node->colorIndex_[1] = newIndex;
						newIndex++;
					}
				}
			}
			else { //node->trackedIndex_==0
				node->colorIndex_[0] = 0;
				node->colorIndex_[1] = 0;
			}
		}
	}
}

void Tracking::setIndexToTrackedStructures() const
{
	size_t newIndex = 1;
	size_t nbSuccessiveNodes;
	Node *node, *oldNode;
	for (size_t time = 0; time < nbTimes_; time++) {
		for (size_t indexNode = 0; indexNode < timeNodes_[time]->nNodes_; indexNode++) {
			if (timeNodes_[time]->nodes_[indexNode]->past_ == nullptr) {
				for (node = timeNodes_[time]->nodes_[indexNode], oldNode = nullptr, nbSuccessiveNodes = 0;
				     node != nullptr && nbSuccessiveNodes < parameters_->minLife_; node = node->future_, nbSuccessiveNodes++) {
					oldNode = node;
					node->trackedIndex_ = 0;
				}
				bool isNewIndex = false;
				if (nbSuccessiveNodes >= parameters_->minLife_) { isNewIndex = true; } //it was == 
				if (timeNodes_[time]->nodes_[indexNode]->fissionPast_ != nullptr) { isNewIndex = true; }
				if (timeNodes_[time]->nodes_[indexNode]->fusion1Past_ != nullptr) { isNewIndex = true; }
				if (timeNodes_[time]->nodes_[indexNode]->blinkPast_ != nullptr) { isNewIndex = true; }
				if (oldNode != nullptr) {
					if (oldNode->fission1Future_ != nullptr) { isNewIndex = true; }
					if (oldNode->fusionFuture_ != nullptr) { isNewIndex = true; }
					if (oldNode->blinkFuture_ != nullptr) { isNewIndex = true; }
				}

				if (isNewIndex) {
					timeNodes_[time]->nodes_[indexNode]->trackedIndex_ = newIndex;
					newIndex++;
				}
			}
			else { timeNodes_[time]->nodes_[indexNode]->trackedIndex_ = timeNodes_[time]->nodes_[indexNode]->past_->trackedIndex_; }
		}
	}
}

/********************************************************************************************************************
int giveResult() just call the function to display the result
********************************************************************************************************************/
void Tracking::GiveResult() const
{
	setIndexToTrackedStructures();
	setColorToTrackedStructures();
	printConnectionsASCII_trc();
}

/********************************************************************************************************************
computeConstantes() compute by convergence the values
parameters_->stdIntensity
parameters_->stdMean
parameters_->stdDistance
********************************************************************************************************************/
void Tracking::computeConstantes() const
{
	parameters_->stdIntensity_ = 1000;
	parameters_->stdMean_ = 1000;
	parameters_->stdDistance_ = 1000;
	double stdMeanComputing = 0;
	double stdDistanceComputing = 0;
	double stdIntensityComputing = 0;


	for (size_t k = 0; k < 20; ++k) {
		long countStdInt = 0, countStdMean = 0, coutStdDist = 0;
		stdMeanComputing = 0;
		stdDistanceComputing = 0;
		stdIntensityComputing = 0;


		for (size_t t = 0; t < nbTimes_ - 1; ++t) {
			timeNodes_[t]->ComputeConstantSigmas(stdIntensityComputing, countStdInt, stdMeanComputing, countStdMean, stdDistanceComputing, coutStdDist, false);
		}

		if (stdIntensityComputing > 1e-8 && countStdInt > 0) { parameters_->stdIntensity_ = stdIntensityComputing / countStdInt; } //*1.0275;//3.4352;
		parameters_->stdIntensity_ = sqrt(parameters_->stdIntensity_);

		if (stdMeanComputing > 1e-8 && countStdMean > 0) { parameters_->stdMean_ = stdMeanComputing / countStdMean; } //*1.0275;//3.4352;
		parameters_->stdMean_ = sqrt(parameters_->stdMean_);

		if (stdDistanceComputing > 1e-8 && coutStdDist > 0) { parameters_->stdDistance_ = stdDistanceComputing / coutStdDist; } //*1.0275;//3.4352;
		parameters_->stdDistance_ = sqrt(parameters_->stdDistance_);
	}
	parameters_->stdDistance_ /= 3; //attention au /3 c'est just pour augmenter l'influence de la distance
}

/********************************************************************************************************************
updateConstantes() ones the Tracking in compute, compute the constant
********************************************************************************************************************/
void Tracking::updateConstantes() const
{
	long cpt = 0;
	double covID = 0;
	double covDM = 0;


	parameters_->stdDistance_ = 0;
	parameters_->stdIntensity_ = 0;
	parameters_->stdMean_ = 0;
	parameters_->covMeanIntensity_ = 0;
	for (size_t t = 0; t < nbTimes_; ++t) {
		const TimeNodes* local = timeNodes_[t];
		for (size_t n = 0; n < local->nNodes_; n++) {
			const Node* node = local->nodes_[n];
			if (node->future_ != nullptr) {
				parameters_->stdDistance_ += node->Distance(node->future_) * node->Distance(node->future_);
				parameters_->stdIntensity_ += (node->info_.intensity_ - node->future_->info_.intensity_)
					* (node->info_.intensity_ - node->future_->info_.intensity_);
				parameters_->stdMean_ += (node->info_.mean_ - node->future_->info_.mean_)
					* (node->info_.mean_ - node->future_->info_.mean_);
				covID += (node->info_.intensity_ - node->future_->info_.intensity_) * node->Distance(node->future_);
				parameters_->covMeanIntensity_ += (node->info_.intensity_ - node->future_->info_.intensity_)
					* (node->info_.mean_ - node->future_->info_.mean_);
				covDM += node->Distance(node->future_) * (node->info_.mean_ - node->future_->info_.mean_);
				cpt++;
			}
		}
	}
	parameters_->stdDistance_ /= cpt;
	parameters_->stdIntensity_ /= cpt;
	parameters_->stdMean_ /= cpt;
	parameters_->covMeanIntensity_ /= cpt;
	parameters_->stdIntensity_ = sqrt(parameters_->stdIntensity_);
	parameters_->stdDistance_ = sqrt(parameters_->stdDistance_);
	parameters_->stdMean_ = sqrt(parameters_->stdMean_);
}

/********************************************************************************************************************
getLikenessGlobal() compute the totalLikeness of hte system
********************************************************************************************************************/
double Tracking::getLikenessGlobal() const
{
	double res = 0.;
	for (size_t t = 0; t < nbTimes_; ++t) { res += timeNodes_[t]->GetLikenessTime(nbTimes_); }
	return res;
}
