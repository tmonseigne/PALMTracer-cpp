#include "Tracking.hpp"
#include <algorithm>
#include <cmath>

#include "Constants.hpp"

namespace palm_tracer::tracking {
//-------------------------------------------------------------------------------------------------
size_t Tracking::Processing()
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

		if (step == 0) {
			computeConstants();
			for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].MakeAllInitializationConnections(); }
		}
		else { updateConstants(); }

		size_t cpt = 0;
		while (cpt < 10) {
			if (step == Parameters::N_STEPS - 1) { error *= 1. - decrease_; }

			double likeness = 0.0;
			for (size_t i = 0; i < Parameters::N_ITERATIONS; ++i) {
				for (size_t t = 0; t < nbTimes_ - 1; ++t) { likeness += timeNodes_[t].TryRandomChange(error); }
			}

			if (std::abs(likeness) < 1e-8) { ++cpt; }
			else { cpt = 0; }
		}
	}

	// Résultat
	setIndexToTrackedStructures();
	setColorToTrackedStructures();
	return printConnectionsASCII_trc();
}

//-------------------------------------------------------------------------------------------------
size_t Tracking::printConnectionsASCII_trc()
{
	size_t idx = 0;
	for (const auto& T : timeNodes_) {
		for (const auto& N : T.nodes_) {
			if (N.past_ == nullptr && N.trackedIndex_ > 0) {
				for (const Node* node = &N; node != nullptr; node = node->future_) {
					tracks_[idx++] = static_cast<double>(N.trackedIndex_); // Track
					tracks_[idx++] = static_cast<double>(node->time_ + 1); // Plane
					tracks_[idx++] = node->id_;                            // Id
					tracks_[idx++] = node->x_;                             // X
					tracks_[idx++] = node->y_;                             // Y
					tracks_[idx++] = node->z_;                             // Z
					tracks_[idx++] = node->intensity_;                     // Intensité
					tracks_[idx++] = node->surf_;                          // Surface

					// Nombre de trajectoires
					nTraj_ = N.trackedIndex_;
				}
			}
		}
	}
	return idx;
}

//-------------------------------------------------------------------------------------------------
void Tracking::readPointsList()
{
	timeNodes_.reserve(nbTimes_);
	params_->maxX_ = 0;
	params_->maxY_ = 0;
	params_->maxZ_ = 0;

	for (size_t t = 0, index = 0; t < nbTimes_; ++t) {
		// Calcul du nombre de noeuds par frame
		size_t nNodes = 0;
		for (size_t i = index; ; i += N_IN_PARAMS) {
			if (static_cast<int>(points_[i + 5]) != -1 && static_cast<int>(points_[i + 5]) != 0) { nNodes++; } // surface
			else { break; }
		}

		timeNodes_.emplace_back(t, nNodes, params_);
		for (size_t i = 0; i < timeNodes_[t].nNodes_; ++i) {
			const size_t Idx = index + i * N_IN_PARAMS;
			// Gestion du cas "2D"
			Node& str      = timeNodes_[t][i];
			str.id_        = points_[Idx + 0];
			str.x_         = points_[Idx + 1];
			str.y_         = points_[Idx + 2];
			str.z_         = points_[Idx + 3];
			str.intensity_ = points_[Idx + 4];
			str.surf_      = points_[Idx + 5];
			str.mean_      = points_[Idx + 4] / points_[Idx + 5];

			params_->maxX_ = std::max(params_->maxX_, str.x_);
			params_->maxY_ = std::max(params_->maxY_, str.y_);
		}

		// Tri des nœuds par intensité décroissante pour un instant t.
		std::ranges::sort(timeNodes_[t].nodes_, [](const Node& a, const Node& b) { return a.intensity_ > b.intensity_; });
		index += (nNodes + 1) * N_IN_PARAMS;
	}

	// Pré-calcul de la distance aux bords
	for (auto& t : timeNodes_) { for (auto& node : t.nodes_) { node.ComputeBorderDistance(); } }
}

//-------------------------------------------------------------------------------------------------
void Tracking::setColorToTrackedStructures()
{
	int idx = 1;
	for (auto& t : timeNodes_) {
		for (auto& node : t.nodes_) {
			if (node.trackedIndex_ != 0) {
				if (node.past_ == nullptr) {
					node.colorIndex_[0] = idx;
					node.colorIndex_[1] = idx;
					idx++;
				}
				else { // node.past_ != nullptr
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

//-------------------------------------------------------------------------------------------------
void Tracking::setIndexToTrackedStructures()
{
	size_t idx = 1, nNodes;
	Node* node;
	for (auto& t : timeNodes_) {
		for (size_t i = 0; i < t.nNodes_; ++i) {
			if (t.nodes_[i].past_ == nullptr) {
				for (node = &t.nodes_[i], nNodes = 0; node != nullptr && nNodes < minLife_; node = node->future_, ++nNodes) { node->trackedIndex_ = 0; }
				if (nNodes >= minLife_) {
					t.nodes_[i].trackedIndex_ = idx;
					idx++;
				}
			}
			else { t.nodes_[i].trackedIndex_ = t.nodes_[i].past_->trackedIndex_; }
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Tracking::computeConstants()
{
	params_->stdIntensity_ = 1000;
	params_->stdMean_      = 1000;
	params_->stdDistance_  = 1000;

	for (size_t k = 0; k < 20; ++k) {
		size_t nInt = 0, nMean    = 0, nDist     = 0;
		double mean = 0, distance = 0, intensity = 0;

		for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].ComputeConstantSigmas(intensity, nInt, mean, nMean, distance, nDist); }

		if (intensity > 1e-8 && nInt > 0) { params_->stdIntensity_ = intensity / static_cast<double>(nInt); }
		params_->stdIntensity_ = std::sqrt(params_->stdIntensity_);

		if (mean > 1e-8 && nMean > 0) { params_->stdMean_ = mean / static_cast<double>(nMean); }
		params_->stdMean_ = std::sqrt(params_->stdMean_);

		if (distance > 1e-8 && nDist > 0) { params_->stdDistance_ = distance / static_cast<double>(nDist); }
		params_->stdDistance_ = std::sqrt(params_->stdDistance_);
	}
	params_->stdDistance_ /= 3; // Augmente l'influence de la distance (pondération empirique).
	for (auto& t : timeNodes_) { for (auto& n : t.nodes_) { n.UpdateConstants(); } }
}

//-------------------------------------------------------------------------------------------------
void Tracking::updateConstants()
{
	size_t cpt                 = 0;
	params_->stdDistance_      = 0;
	params_->stdIntensity_     = 0;
	params_->stdMean_          = 0;
	params_->covMeanIntensity_ = 0;

	for (const auto& T : timeNodes_) {
		for (const auto& Node : T.nodes_) {
			if (Node.future_ != nullptr) {
				params_->stdDistance_ += Node.Distance(*Node.future_) * Node.Distance(*Node.future_);
				params_->stdIntensity_ += (Node.intensity_ - Node.future_->intensity_) * (Node.intensity_ - Node.future_->intensity_);
				params_->stdMean_ += (Node.mean_ - Node.future_->mean_) * (Node.mean_ - Node.future_->mean_);
				params_->covMeanIntensity_ += (Node.intensity_ - Node.future_->intensity_) * (Node.mean_ - Node.future_->mean_);
				cpt++;
			}
		}
	}

	// Sécurise une éventuelle division par zéro si aucune liaison n'a été trouvée.
	if (cpt == 0) {
		params_->stdIntensity_     = 0;
		params_->stdDistance_      = 0;
		params_->stdMean_          = 0;
		params_->covMeanIntensity_ = 0;
		return;
	}

	params_->stdDistance_ /= static_cast<double>(cpt);
	params_->stdIntensity_ /= static_cast<double>(cpt);
	params_->stdMean_ /= static_cast<double>(cpt);
	params_->covMeanIntensity_ /= static_cast<double>(cpt);
	params_->stdIntensity_ = std::sqrt(params_->stdIntensity_);
	params_->stdDistance_  = std::sqrt(params_->stdDistance_);
	params_->stdMean_      = std::sqrt(params_->stdMean_);
	for (auto& t : timeNodes_) { for (auto& n : t.nodes_) { n.UpdateConstants(); } }
}
} // namespace palm_tracer::tracking
