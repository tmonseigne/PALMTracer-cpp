///-------------------------------------------------------------------------------------------------
///
/// \file TimeNodes.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/
#pragma once

#include "Node.hpp"
#include "Parameters.hpp"

namespace palm_tracer {
class TimeNodes {
protected:
	size_t time_;

public:
	Parameters* params_;
	std::vector<Node> nodes_;
	size_t nNodes_;
	TimeNodes(size_t time, size_t nNodes, Parameters* params);
	virtual ~TimeNodes() { nodes_.clear(); }

	Node& operator[](const size_t i) { return nodes_[i]; }
	void LinkNeighborsBeforeAndAfter(TimeNodes& futures);
	void SortByIntensityAllNeighbors();

	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const;
	void MakeAllInitializationConnections() { for (auto& node : nodes_) { node.MakeInitializationConnection(); } }
	double TryRandomChange(const double limit) { return nNodes_ == 0 ? 0 : nodes_[params_->RandomModulo(nNodes_)].MakeChange(limit); }
};
} // namespace palm_tracer
