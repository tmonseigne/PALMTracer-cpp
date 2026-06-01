///-------------------------------------------------------------------------------------------------
///
/// \file TrackingPoint.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 28/05/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include "Point.hpp"

namespace palm_tracer::tracking {
static constexpr size_t INVALID_INDEX = std::numeric_limits<size_t>::max();

/// <summary></summary>
struct TrackingPoint {
	Point point;
	size_t next     = INVALID_INDEX;
	size_t previous = INVALID_INDEX;
	size_t trackId  = 0;

	/// <summary></summary>
	/// <param name="other"></param>
	/// <returns></returns>
	[[nodiscard]] double GetSqDistance(const TrackingPoint& other) const
	{
		const auto P = point.posPx_ - other.point.posPx_;
		return P.x * P.x + P.y * P.y;
	}
};

/// <summary></summary>
struct PlaneRange {
	size_t begin = 0;
	size_t end   = 0;
};

/// <summary></summary>
struct Candidate {
	size_t from      = INVALID_INDEX;
	size_t to        = INVALID_INDEX;
	double distance2 = 0.0;
};

/// <summary></summary>
struct Track {
	size_t id = 0;
	std::vector<Point> points; ///< Liste des points de la trajectoire.
};
} // namespace palm_tracer::tracking
