///-------------------------------------------------------------------------------------------------
///
/// \file Point.hpp
/// \brief Définition de la classe \c Point représentant un élément d'une trajectoire.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 30/09/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include "Structs.hpp"


namespace palm_tracer {
/// <summary> Représente un point élémentaire d'une trajectoire.
/// Un point regroupe les informations spatiales (plan, coordonnées), ainsi que des caractéristiques photométriques (intensité, surface).</summary>
class Point {
public:
	size_t id_        = 0;                        ///< Identifiant du point dans le tableau initial ou source.
	size_t plane_     = 0;                        ///< Numéro de plan (frame) associé au point.
	Vec3D posPx_      = {.x = 0, .y = 0, .z = 0}; ///< Position 3D du point (coordonnées en pixels).
	Vec3D pos_        = {.x = 0, .y = 0, .z = 0}; ///< Position 3D du point (coordonnées en unités réelles).
	double intensity_ = 0,                        ///< Itensité totale mesurée dans la région d'intérêt (somme des valeurs de pixels).
	       surface_   = 0;                        ///< Surface de la région associée au point (aire en pixels²).

	Point()  = default; ///< Constructeur par défaut de la classe <c>Point</c>. Initialise tous les attributs à zéro.
	~Point() = default; ///< Destructeur par défaut de la classe <c>Point</c>.

	/// <summary>Constructeur paramétré de la classe <see cref="Point"/>.</summary>
	/// <param name="id">Identifiant unique du point dans le tableau initial.</param>
	/// <param name="plane">Numéro du plan (frame) auquel le point appartient.</param>
	/// <param name="pPx">Position spatiale (coordonnées X, Y, Z en pixel).</param>
	/// <param name="p">Position spatiale (coordonnées X, Y, Z).</param>
	/// <param name="intensity">Intensité totale de la région associée.</param>
	/// <param name="surface">Surface de la région associée.</param>
	explicit Point(const size_t id, const size_t plane, const Vec3D& pPx, const Vec3D& p, const double intensity, const double surface)
		: id_(id), plane_(plane), posPx_(pPx), pos_(p), intensity_(intensity), surface_(surface) { }

	/// <summary>Calcule la distance euclidienne entre ce point et un autre.</summary>
	/// <param name="other">Autre point avec lequel calculer la distance.</param>
	/// <returns>Distance euclidienne (même unité que les coordonnées de <c>pos_</c>).</returns>
	[[nodiscard]] double GetDistance(const Point& other) const { return (pos_ - other.pos_).Norm(); }

	/// <summary>Surcharge de l’opérateur de flux de sortie pour afficher un point de trajectoire.</summary>
	/// <param name="os">Flux de sortie.</param>
	/// <param name="p">Objet <c>point</c> à afficher.</param>
	/// <returns>Référence vers le flux de sortie, pour chaînage.</returns>
	friend std::ostream& operator<<(std::ostream& os, const Point& p)
	{
		os << "[" << p.id_ << ", " << p.plane_ << ", " << p.posPx_ << ", " << p.intensity_ << ", " << p.surface_ << "]";
		return os;
	}
};
} // namespace palm_tracer
