///-------------------------------------------------------------------------------------------------
///
/// \file Point.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include "Structs.hpp"


namespace palm_tracer {
/// <summary>Classe manipulant un point d'une trajectoire.</summary>
class Point {
public:
	size_t id_        = 0;                        ///< Index du point dans le tableau initial.
	size_t plane_     = 0;                        ///< Plan du point.
	Vec3D pos_        = {.x = 0, .y = 0, .z = 0}; ///< Position.
	double intensity_ = 0,                        ///< Intensité totale dans la région.
	       surface_   = 0;                        ///< Surface de la région.

	Point()  = default; ///< Constructeur de la classe <c>Point</c>.
	~Point() = default; ///< Destructeur de la classe <c>Point</c>.

	/// <summary>Constructeur de la classe <see cref="Point"/> avec paramètres.</summary>
	/// <param name="id">Index du point dans le tableau initial.</param>
	/// <param name="plane">Plan du point.</param>
	/// <param name="p">Position.</param>
	/// <param name="intensity">Intensité totale dans la région.</param>
	/// <param name="surface">Surface de la région.</param>
	explicit Point(const size_t id, const size_t plane, const Vec3D& p, const double intensity, const double surface)
		: id_(id), plane_(plane), pos_(p), intensity_(intensity), surface_(surface) { }

	/// <summary>Calcule la distance euclidienne entre le point et un autre.</summary>
	/// <param name="other">The other.</param>
	/// <returns></returns>
	[[nodiscard]] double GetDistance(const Point& other) const { return (pos_ - other.pos_).Norm(); }
};
} // namespace palm_tracer
