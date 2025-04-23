/*******************************************************************************************************************
									Tracking Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/

/********************************************************************************************************************
Class Tracking contains several Tracking objets for each time (nbTimes).
********************************************************************************************************************/


#include "Tracking.h"
#include <algorithm>

void Tracking::Processing()
{
	// Initialisation
	readPointsList();
	linkNeighborsBeforeAndAfter();
	sortByIntensityAllNeighbors();
	if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1) { linkNeighborsBlinkBeforeAndAfter(); }

	// Process
	double error;
	constexpr size_t MaxStep = 5;

	for (size_t step = 0; step < MaxStep; ++step) {
		if (step == MaxStep - 1) { error = 20; }
		else { error = 0; }
		long cpt = 0;

		if (step == 0) {
			computeConstantes();
			for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].MakeAllInitializationConnections(); }
		}
		else { updateConstantes(); }

		int nbInt = 0;
		int bcl = 0;
		do {
			if (step == MaxStep - 1) { error *= 1. - decrease_; }

			nbInt++;

			double likeness = 0;
			for (size_t i = 0; i < Parameters::N_ITERATIONS; ++i) {
				for (size_t t = 0; t < nbTimes_ - 1; ++t) { likeness += timeNodes_[t].TryRandomChange(error); }
			}

			if (fabs(likeness) < 1e-8) { cpt++; }
			else { cpt = 0; }

			bcl++;
		}
		while (cpt < 10);
	}

	// Result
	setIndexToTrackedStructures();
	setColorToTrackedStructures();
	printConnectionsASCII_trc();
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
void Tracking::printConnectionsASCII_trc()
{
	size_t idx = 0;
	for (const auto& T : timeNodes_) {
		for (const auto& N : T.nodes_) {
			if (N.past_ == nullptr && N.trackedIndex_ > 0) {
				for (const Node* nodeConnection = &N; nodeConnection != nullptr; nodeConnection = nodeConnection->future_) {
					tracks_[idx++] = static_cast<double>(N.trackedIndex_);
					tracks_[idx++] = static_cast<double>(nodeConnection->time_ + 1);
					tracks_[idx++] = nodeConnection->info_.x_;
					tracks_[idx++] = nodeConnection->info_.y_;
					tracks_[idx++] = nodeConnection->info_.surf_;
					tracks_[idx++] = nodeConnection->info_.intensity_;
					tracks_[idx++] = nodeConnection->info_.z_;
					tracks_[idx++] = nodeConnection->info_.color_;
					tracks_[idx++] = nodeConnection->info_.molId_;

					// The number of trajectory
					nTraj_ = N.trackedIndex_;
				}
			}
		}
	}
}

/********************************************************************************************************************
readPointsList() reads the 'pointsList' specified in params_ and allocates the corresponding TimeNodes
********************************************************************************************************************/
void Tracking::readPointsList()
{
	timeNodes_.reserve(nbTimes_);
	params_->maxX_ = 0;
	params_->maxY_ = 0;
	params_->maxZ_ = 0;

	for (size_t t = 0, index = 0; t < nbTimes_; ++t) {
		// Computing the number of nodes by frame
		size_t nNodes = 0;
		for (size_t i = index; ; i += Parameters::N_SEG_PARAMS) {
			if (static_cast<int>(points_[i + 9]) != -1 && static_cast<int>(points_[i + 9]) != 0) { nNodes++; } // surface
			else { break; }
		}

		timeNodes_.emplace_back(t, nNodes, params_);
		for (size_t i = 0; i < timeNodes_[t].nNodes_; ++i) {
			const size_t Idx = index + i * Parameters::N_SEG_PARAMS;
			// Just for the case of '2D'
			InfoNode& str = timeNodes_[t][i].info_;
			// version correcte
			str.x_ = points_[Idx + 3];
			str.y_ = points_[Idx + 4];
			str.z_ = points_[Idx + 10];

			str.color_ = points_[Idx + 11];
			str.molId_ = points_[Idx + 12];

			str.intensity_ = points_[Idx + 8];

			str.surf_ = points_[Idx + 9];
			str.mean_ = points_[Idx + 8] / points_[Idx + 9];

			params_->maxX_ = std::max(params_->maxX_, timeNodes_[t][i].info_.x_);
			params_->maxY_ = std::max(params_->maxY_, timeNodes_[t][i].info_.y_);
		}

		//Sorting the nodes by the intensity for one time
		size_t indexMax = 0;
		for (size_t i = 0; i < timeNodes_[t].nNodes_ - 1; ++i) {
			double maxIntensity = -1.;
			for (size_t j = i; j < timeNodes_[t].nNodes_; ++j) {
				if (maxIntensity < timeNodes_[t][j].info_.intensity_) {
					maxIntensity = timeNodes_[t][j].info_.intensity_;
					indexMax = j;
				}
			}
			std::swap(timeNodes_[t].nodes_[indexMax], timeNodes_[t].nodes_[i]);
		}

		index += (nNodes + 1) * Parameters::N_SEG_PARAMS; //(A vérifier)
	}
}

/********************************************************************************************************************
void linkNeighborsBlinkBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of 
all times for blinking purpose
********************************************************************************************************************/
void Tracking::linkNeighborsBlinkBeforeAndAfter()
{
	for (size_t t = 0; t < nbTimes_; ++t) {
		for (size_t dt = 2; dt <= Parameters::N_TEMPORAL_WINDOWS; ++dt) {
			if (t + dt < nbTimes_) { timeNodes_[t].LinkNeighborsBlinkBeforeAndAfter(timeNodes_[t + dt]); }
		}
	}
}

/********************************************************************************************************************
void setColorToTrackedStructures() set a color index to each tracked structures
********************************************************************************************************************/
void Tracking::setColorToTrackedStructures()
{
	Node *nodePastFission, *nodePastFission2;
	int idx = 1;
	for (auto& T : timeNodes_) {
		for (auto& node : T.nodes_) {
			if (node.trackedIndex_ != 0) {
				if (node.past_ == nullptr) {
					if (node.blinkPast_ == nullptr) {
						if (node.fusion1Past_ != nullptr) {
							if (node.Likeness(node.fusion1Past_) > node.Likeness(node.fusion2Past_)) {
								node.colorIndex_[0] = node.fusion1Past_->colorIndex_[0];
								node.colorIndex_[1] = node.fusion2Past_->colorIndex_[0];
							}
							else {
								node.colorIndex_[0] = node.fusion2Past_->colorIndex_[0];
								node.colorIndex_[1] = node.fusion1Past_->colorIndex_[0];
							}
						}
						else { //node.fusion1Past_==NULL
							if (node.fissionPast_ != nullptr) {
								Node* node2 = &node == node.fissionPast_->fission1Future_
									              ? node.fissionPast_->fission2Future_
									              : node.fissionPast_->fission1Future_;

								nodePastFission2 = node.fissionPast_;
								for (nodePastFission = node.fissionPast_; nodePastFission != nullptr; nodePastFission = nodePastFission->past_) {
									nodePastFission2 = nodePastFission;
								}
								if (nodePastFission2->fusion1Past_ != nullptr) {
									const Node* nodePastFission3 = nodePastFission2->fusion1Past_;
									const Node* nodePastFission4 = nodePastFission2->fusion2Past_;
									if (node.Likeness(nodePastFission3) + node2->Likeness(nodePastFission4) > node2->Likeness(nodePastFission3) + node.
									    Likeness(nodePastFission4)) {
										node.colorIndex_[0] = nodePastFission3->colorIndex_[0];
										node.colorIndex_[1] = nodePastFission3->colorIndex_[0];
										node2->colorIndex_[0] = nodePastFission4->colorIndex_[0];
										node2->colorIndex_[1] = nodePastFission4->colorIndex_[0];
									}
									else {
										node2->colorIndex_[0] = nodePastFission3->colorIndex_[0];
										node2->colorIndex_[1] = nodePastFission3->colorIndex_[0];
										node.colorIndex_[0] = nodePastFission4->colorIndex_[0];
										node.colorIndex_[1] = nodePastFission4->colorIndex_[0];
									}
								}
								else {
									if (node.fissionPast_->Likeness(&node) > node.fissionPast_->Likeness(node2)) {
										node.colorIndex_[0] = node.fissionPast_->colorIndex_[0];
										node.colorIndex_[1] = node.fissionPast_->colorIndex_[0];
									}
									else {
										node.colorIndex_[0] = node.fissionPast_->colorIndex_[1];
										node.colorIndex_[1] = node.fissionPast_->colorIndex_[1];
									}
								}
							}
							else { //node.fissionPast_==NULL  and node.fusion1Past_==NULL
								if (node.fission1Future_ != nullptr) {
									node.colorIndex_[0] = idx;
									idx++;
									node.colorIndex_[1] = idx;
									idx++;
								}
								else {
									node.colorIndex_[0] = idx;
									node.colorIndex_[1] = idx;
									idx++;
								}
							}
						}
					}
					else {
						node.colorIndex_[0] = node.blinkPast_->colorIndex_[0];
						node.colorIndex_[1] = node.blinkPast_->colorIndex_[0];
					}
				}
				else { //node.past_!=NULL
					node.colorIndex_[0] = node.past_->colorIndex_[0];
					node.colorIndex_[1] = node.past_->colorIndex_[0];
				}
				if (node.fission1Future_ != nullptr) {
					nodePastFission2 = &node;
					for (nodePastFission = &node; nodePastFission != nullptr; nodePastFission = nodePastFission->past_) { nodePastFission2 = nodePastFission; }
					if (nodePastFission2->fusion1Past_ != nullptr) {
						node.colorIndex_[0] = nodePastFission2->fusion1Past_->colorIndex_[0];
						node.colorIndex_[1] = nodePastFission2->fusion2Past_->colorIndex_[0];
					}
					else {
						node.colorIndex_[1] = idx;
						idx++;
					}
				}
			}
			else { //node.trackedIndex_==0
				node.colorIndex_[0] = 0;
				node.colorIndex_[1] = 0;
			}
		}
	}
}

void Tracking::setIndexToTrackedStructures()
{
	size_t idx = 1, nNodes;
	Node *node, *oldNode;
	for (auto& T : timeNodes_) {
		for (size_t i = 0; i < T.nNodes_; ++i) {
			if (T.nodes_[i].past_ == nullptr) {
				for (node = &T.nodes_[i], oldNode = nullptr, nNodes = 0;
				     node != nullptr && nNodes < minLife_; node = node->future_, ++nNodes) {
					oldNode = node;
					node->trackedIndex_ = 0;
				}
				bool change = false;
				if (nNodes >= minLife_) { change = true; } //it was == 
				if (T.nodes_[i].fissionPast_ != nullptr) { change = true; }
				if (T.nodes_[i].fusion1Past_ != nullptr) { change = true; }
				if (T.nodes_[i].blinkPast_ != nullptr) { change = true; }
				if (oldNode != nullptr) {
					if (oldNode->fission1Future_ != nullptr) { change = true; }
					if (oldNode->fusionFuture_ != nullptr) { change = true; }
					if (oldNode->blinkFuture_ != nullptr) { change = true; }
				}

				if (change) {
					T.nodes_[i].trackedIndex_ = idx;
					idx++;
				}
			}
			else { T.nodes_[i].trackedIndex_ = T.nodes_[i].past_->trackedIndex_; }
		}
	}
}

/********************************************************************************************************************
computeConstantes() compute by convergence the values
params_->stdIntensity
params_->stdMean
params_->stdDistance
********************************************************************************************************************/
void Tracking::computeConstantes() const
{
	params_->stdIntensity_ = 1000;
	params_->stdMean_ = 1000;
	params_->stdDistance_ = 1000;

	for (size_t k = 0; k < 20; ++k) {
		size_t nInt = 0, nMean = 0, nDist = 0;
		double mean = 0, distance = 0, intensity = 0;

		for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].ComputeConstantSigmas(intensity, nInt, mean, nMean, distance, nDist, false); }

		if (intensity > 1e-8 && nInt > 0) { params_->stdIntensity_ = intensity / static_cast<double>(nInt); } //*1.0275;//3.4352;
		params_->stdIntensity_ = sqrt(params_->stdIntensity_);

		if (mean > 1e-8 && nMean > 0) { params_->stdMean_ = mean / static_cast<double>(nMean); } //*1.0275;//3.4352;
		params_->stdMean_ = sqrt(params_->stdMean_);

		if (distance > 1e-8 && nDist > 0) { params_->stdDistance_ = distance / static_cast<double>(nDist); } //*1.0275;//3.4352;
		params_->stdDistance_ = sqrt(params_->stdDistance_);
	}
	params_->stdDistance_ /= 3; //attention au /3 c'est juste pour augmenter l'influence de la distance
}

/********************************************************************************************************************
updateConstantes() ones the Tracking in compute, compute the constant
********************************************************************************************************************/
void Tracking::updateConstantes() const
{
	long cpt = 0;
	params_->stdDistance_ = 0;
	params_->stdIntensity_ = 0;
	params_->stdMean_ = 0;
	params_->covMeanIntensity_ = 0;

	for (const auto& T : timeNodes_) {
		for (const auto& Node : T.nodes_) {
			if (Node.future_ != nullptr) {
				params_->stdDistance_ += Node.Distance(*Node.future_) * Node.Distance(*Node.future_);
				params_->stdIntensity_ += (Node.info_.intensity_ - Node.future_->info_.intensity_)
					* (Node.info_.intensity_ - Node.future_->info_.intensity_);
				params_->stdMean_ += (Node.info_.mean_ - Node.future_->info_.mean_)
					* (Node.info_.mean_ - Node.future_->info_.mean_);
				params_->covMeanIntensity_ += (Node.info_.intensity_ - Node.future_->info_.intensity_)
					* (Node.info_.mean_ - Node.future_->info_.mean_);
				cpt++;
			}
		}
	}
	params_->stdDistance_ /= cpt;
	params_->stdIntensity_ /= cpt;
	params_->stdMean_ /= cpt;
	params_->covMeanIntensity_ /= cpt;
	params_->stdIntensity_ = sqrt(params_->stdIntensity_);
	params_->stdDistance_ = sqrt(params_->stdDistance_);
	params_->stdMean_ = sqrt(params_->stdMean_);
}
