///-------------------------------------------------------------------------------------------------
///
/// \file Structs.hpp
/// \brief  Fichier de structures necessaire aux algorithmes.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

namespace palm_tracer {
/// <summary> Structure représentant un pixel avec ses coordonnées (x, y) et une valeur scalaire associée (v).</summary>
struct Pixel {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y; ///< Coordonnée entière verticale du pixel (axe Y).
	double v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).

	// Pour un tri par valeur
	bool operator<(const Pixel& other) const { return v > other.v; }
};

/// <summary> Structure représentant un pixel avec ses coordonnées (x, y) et une valeur scalaire associée (v).</summary>
struct Voxel {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y, ///< Coordonnée entière verticale du pixel (axe Y).
	       z; ///< Coordonnée entière verticale du pixel (axe Y).
	double v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};

/// <summary> Structure représentant un pixel avec ses coordonnées (x, y) en flottant et une valeur scalaire associée (v).</summary>
struct Point2D {
	double x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y, ///< Coordonnée entière verticale du pixel (axe Y).
	       v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};

/// <summary> Structure représentant un pixel avec ses coordonnées (x, y) en flottant et une valeur scalaire associée (v).</summary>
struct Point3D {
	double x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y, ///< Coordonnée entière verticale du pixel (axe Y).
	       z, ///< Coordonnée entière verticale du pixel (axe Y).
	       v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};
} // namespace palm_tracer
