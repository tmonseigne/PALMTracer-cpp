#include "TimeNodes.hpp"
#include <algorithm>
#include <ranges>

namespace palm_tracer::tracking {
//-------------------------------------------------------------------------------------------------
TimeNodes::TimeNodes(const size_t time, const size_t nNodes, Parameters* params)
	: time_(time), params_(params), nNodes_(nNodes)
{
	nodes_.reserve(nNodes_);
	for (size_t i = 0; i < nNodes_; ++i) { nodes_.emplace_back(time_, params_); }
}

//-------------------------------------------------------------------------------------------------
void TimeNodes::LinkNeighborsBeforeAndAfter(TimeNodes& futures)
{
	for (auto& present : nodes_) {
		for (auto& future : futures.nodes_) {
			if (present.Distance(future) < params_->maxSpeed_) {
				present.neighbors_[1].emplace_back(&future);
				future.neighbors_[0].emplace_back(&present);
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
void TimeNodes::SortByIntensityAllNeighbors()
{
	for (auto& n : nodes_) {
		for (auto& neightbors : n.neighbors_) { std::ranges::sort(neightbors, [](const Node* a, const Node* b) { return a->intensity_ > b->intensity_; }); }
	}
}


//-------------------------------------------------------------------------------------------------
void TimeNodes::ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const
{
	for (const auto& N : nodes_) { N.ComputeConstantSigmas(intensity, nInt, mean, nMean, distance, nDist); }
}
} // namespace palm_tracer::tracking
