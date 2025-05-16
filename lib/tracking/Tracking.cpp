/*******************************************************************************************************************
									Tracking Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/

/********************************************************************************************************************
Class Tracking contains several Tracking objets for each time (nbTimes).
********************************************************************************************************************/


#include "Tracking.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>

void Tracking::Processing()
{
	// Initialisation
	readPointsList();
	linkNeighborsBeforeAndAfter();
	sortByIntensityAllNeighbors();

	// Process
	double error;

	for (size_t step = 0; step < Parameters::N_STEPS; ++step) {
		if (step == Parameters::N_STEPS - 1) { error = 20; }
		else { error = 0; }
		long cpt = 0;

		if (step == 0) {
			computeConstants();
			for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].MakeAllInitializationConnections(); }
		}
		else { updateConstants(); }

		int nbInt = 0;
		int bcl = 0;
		do {
			if (step == Parameters::N_STEPS - 1) { error *= 1. - decrease_; }

			nbInt++;

			double likeness = 0;
			for (size_t i = 0; i < Parameters::N_ITERATIONS; ++i) {
				for (size_t t = 0; t < nbTimes_ - 1; ++t) { likeness += timeNodes_[t].TryRandomChange(error); }
			}

			if (std::fabs(likeness) < 1e-8) { cpt++; }
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
			if (!N.past_ && N.trackedIndex_ > 0) {
				for (const Node* node = &N; node; node = node->future_) {
					tracks_[idx++] = static_cast<double>(N.trackedIndex_);
					tracks_[idx++] = static_cast<double>(node->time_ + 1);
					tracks_[idx++] = node->x_;
					tracks_[idx++] = node->y_;
					tracks_[idx++] = node->surf_;
					tracks_[idx++] = node->intensity_;
					tracks_[idx++] = node->z_;
					tracks_[idx++] = node->color_;
					tracks_[idx++] = node->molId_;

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
			Node& str = timeNodes_[t][i];
			// version correcte
			str.x_ = points_[Idx + 3];
			str.y_ = points_[Idx + 4];
			str.z_ = points_[Idx + 10];

			str.color_ = points_[Idx + 11];
			str.molId_ = points_[Idx + 12];

			str.intensity_ = points_[Idx + 8];

			str.surf_ = points_[Idx + 9];
			str.mean_ = points_[Idx + 8] / points_[Idx + 9];

			params_->maxX_ = std::max(params_->maxX_, timeNodes_[t][i].x_);
			params_->maxY_ = std::max(params_->maxY_, timeNodes_[t][i].y_);
		}

		//Sorting the nodes by the intensity for one time
		size_t indexMax = 0;
		for (size_t i = 0; i < timeNodes_[t].nNodes_ - 1; ++i) {
			double maxIntensity = -1.;
			for (size_t j = i; j < timeNodes_[t].nNodes_; ++j) {
				if (maxIntensity < timeNodes_[t][j].intensity_) {
					maxIntensity = timeNodes_[t][j].intensity_;
					indexMax = j;
				}
			}
			std::swap(timeNodes_[t].nodes_[indexMax], timeNodes_[t].nodes_[i]);
		}

		index += (nNodes + 1) * Parameters::N_SEG_PARAMS; //(À vérifier)
	}

	//Precompute distance to border
	for (auto& t : timeNodes_) { for (auto& node : t.nodes_) { node.ComputeBorderDistance(); } }
}

/********************************************************************************************************************
void setColorToTrackedStructures() set a color index to each tracked structures
********************************************************************************************************************/
void Tracking::setColorToTrackedStructures()
{
	int idx = 1;
	for (auto& T : timeNodes_) {
		for (auto& node : T.nodes_) {
			if (node.trackedIndex_ != 0) {
				if (!node.past_) {
					node.colorIndex_[0] = idx;
					node.colorIndex_[1] = idx;
					idx++;
				}
				else { //node.past_!=NULL
					node.colorIndex_[0] = node.past_->colorIndex_[0];
					node.colorIndex_[1] = node.past_->colorIndex_[0];
				}
			}
			else {
				node.colorIndex_[0] = 0;
				node.colorIndex_[1] = 0;
			}
		}
	}
}

void Tracking::setIndexToTrackedStructures()
{
	size_t idx = 1, nNodes;
	Node* node;
	for (auto& T : timeNodes_) {
		for (size_t i = 0; i < T.nNodes_; ++i) {
			if (!T.nodes_[i].past_) {
				for (node = &T.nodes_[i], nNodes = 0; node && nNodes < minLife_; node = node->future_, ++nNodes) { node->trackedIndex_ = 0; }
				if (nNodes >= minLife_) {
					T.nodes_[i].trackedIndex_ = idx;
					idx++;
				}
			}
			else { T.nodes_[i].trackedIndex_ = T.nodes_[i].past_->trackedIndex_; }
		}
	}
}

/********************************************************************************************************************
computeConstants() compute by convergence the values
params_->stdIntensity
params_->stdMean
params_->stdDistance
********************************************************************************************************************/
void Tracking::computeConstants()
{
	params_->stdIntensity_ = 1000;
	params_->stdMean_ = 1000;
	params_->stdDistance_ = 1000;

	for (size_t k = 0; k < 20; ++k) {
		size_t nInt = 0, nMean = 0, nDist = 0;
		double mean = 0, distance = 0, intensity = 0;

		for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].ComputeConstantSigmas(intensity, nInt, mean, nMean, distance, nDist); }

		if (intensity > 1e-8 && nInt > 0) { params_->stdIntensity_ = intensity / static_cast<double>(nInt); } //*1.0275;//3.4352;
		params_->stdIntensity_ = sqrt(params_->stdIntensity_);

		if (mean > 1e-8 && nMean > 0) { params_->stdMean_ = mean / static_cast<double>(nMean); } //*1.0275;//3.4352;
		params_->stdMean_ = sqrt(params_->stdMean_);

		if (distance > 1e-8 && nDist > 0) { params_->stdDistance_ = distance / static_cast<double>(nDist); } //*1.0275;//3.4352;
		params_->stdDistance_ = sqrt(params_->stdDistance_);
	}
	params_->stdDistance_ /= 3; //attention au / 3 c'est juste pour augmenter l'influence de la distance
	for (auto& t : timeNodes_) { for (auto& n : t.nodes_) { n.UpdateConstants(); } }
}

/********************************************************************************************************************
updateConstants() ones the Tracking in compute, compute the constant
********************************************************************************************************************/
void Tracking::updateConstants()
{
	long cpt = 0;
	params_->stdDistance_ = 0;
	params_->stdIntensity_ = 0;
	params_->stdMean_ = 0;
	params_->covMeanIntensity_ = 0;

	for (const auto& T : timeNodes_) {
		for (const auto& Node : T.nodes_) {
			if (Node.future_) {
				params_->stdDistance_ += Node.Distance(*Node.future_) * Node.Distance(*Node.future_);
				params_->stdIntensity_ += (Node.intensity_ - Node.future_->intensity_) * (Node.intensity_ - Node.future_->intensity_);
				params_->stdMean_ += (Node.mean_ - Node.future_->mean_) * (Node.mean_ - Node.future_->mean_);
				params_->covMeanIntensity_ += (Node.intensity_ - Node.future_->intensity_) * (Node.mean_ - Node.future_->mean_);
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
	for (auto& t : timeNodes_) { for (auto& n : t.nodes_) { n.UpdateConstants(); } }
}
