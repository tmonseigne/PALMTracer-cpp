///-------------------------------------------------------------------------------------------------
///
/// \file TrackingPoint.hpp
/// <summary>Structures de travail pour les points, plages de plans et associations de trajectoires.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 28/05/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <limits>
#include "Point.hpp"

namespace palm_tracer::tracking {
static constexpr size_t INVALID_INDEX = std::numeric_limits<size_t>::max(); ///< Absence d'indice de point ou de liaison.

/// <summary>Point associé à ses liaisons vers les plans précédent et suivant.</summary>
struct TrackingPoint {
	Point point;                     ///< Mesures et positions de la localisation.
	size_t next     = INVALID_INDEX; ///< Indice du point suivant, ou INVALID_INDEX en l'absence de liaison.
	size_t previous = INVALID_INDEX; ///< Indice du point précédent, ou INVALID_INDEX en l'absence de liaison.
	size_t trackId  = 0;             ///< Identifiant de trajectoire, 0 avant attribution.

	/// <summary>Calcule la distance au carré entre deux positions XY en pixels.</summary>
	/// <param name="other">Point auquel comparer la position.</param>
	/// <returns>Distance XY au carré, en pixels carrés ; Z n'intervient pas.</returns>
	[[nodiscard]] double GetSqDistance(const TrackingPoint& other) const
	{
		const auto P = point.posPx_ - other.point.posPx_;
		return (P.x * P.x) + (P.y * P.y);
	}
};

/// <summary>Intervalle semi-ouvert des points d'un même plan dans un tableau trié.</summary>
struct PlaneRange {
	size_t begin = 0; ///< Premier indice inclus.
	size_t end   = 0; ///< Premier indice exclu.
};

/// <summary>Liaison candidate entre deux points de plans consécutifs.</summary>
struct Candidate {
	size_t from      = INVALID_INDEX; ///< Indice du point source.
	size_t to        = INVALID_INDEX; ///< Indice du point cible.
	double distance2 = 0.0;           ///< Distance XY au carré, en pixels carrés.
};

/// <summary>Trajectoire identifiée regroupant les localisations dans l'ordre temporel.</summary>
struct Track {
	size_t id = 0;             ///< Identifiant de trajectoire.
	std::vector<Point> points; ///< Liste des points de la trajectoire.
};
} // namespace palm_tracer::tracking
