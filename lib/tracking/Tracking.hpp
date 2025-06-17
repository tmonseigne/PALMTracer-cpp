///-------------------------------------------------------------------------------------------------
///
/// \file Tracking.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
/*******************************************************************************************************************
									Tracking Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/
#pragma once

#include "Parameters.hpp"
#include "TimeNodes.hpp"

namespace palm_tracer {
class Tracking {
public:
	size_t nTraj_;

	explicit Tracking(const double* points, double* tracks, const size_t minLife, const double decrease, const size_t ntimes, Parameters* params)
		: nTraj_(0), points_(points), tracks_(tracks), decrease_(decrease), minLife_(minLife), nbTimes_(ntimes), params_(params) { }

	~Tracking() { timeNodes_.clear(); }

	size_t Processing();

protected:
	const double* points_;
	double* tracks_;
	double decrease_;
	size_t minLife_, nbTimes_;
	Parameters* params_;
	std::vector<TimeNodes> timeNodes_;

	void updateConstants();
	void computeConstants();
	void linkNeighborsBeforeAndAfter() { for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].LinkNeighborsBeforeAndAfter(timeNodes_[t + 1]); } }
	void sortByIntensityAllNeighbors() { for (size_t t = 0; t < nbTimes_; ++t) { timeNodes_[t].SortByIntensityAllNeighbors(); } }
	void readPointsList();
	size_t printConnectionsASCII_trc();
	void setColorToTrackedStructures();
	void setIndexToTrackedStructures();
};
} // namespace palm_tracer
