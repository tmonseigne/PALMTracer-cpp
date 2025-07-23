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
#include <iostream>

namespace palm_tracer {
/// <summary>Structure représentant un pixel avec ses coordonnées (x, y) et une valeur scalaire associée (v).</summary>
struct Pixel {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y; ///< Coordonnée entière verticale du pixel (axe Y).
	double v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).

	// Pour un tri par valeur
	bool operator<(const Pixel& other) const { return v > other.v; }
};

/// <summary>Structure représentant un pixel avec ses coordonnées (x, y) et une valeur scalaire associée (v).</summary>
struct Voxel {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y, ///< Coordonnée entière verticale du pixel (axe Y).
	       z; ///< Coordonnée entière de profondeur du pixel (axe Z).
	double v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};

/// <summary>Structure représentant un pixel avec ses coordonnées (x, y) en flottant et une valeur scalaire associée (v).</summary>
struct Point2D {
	double x, ///< Coordonnée flottante horizontale du pixel (axe X).
	       y, ///< Coordonnée flottante verticale du pixel (axe Y).
	       v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};

/// <summary>Structure représentant un pixel avec ses coordonnées (x, y) en flottant et une valeur scalaire associée (v).</summary>
struct Point3D {
	double x, ///< Coordonnée flottante horizontale du pixel (axe X).
	       y, ///< Coordonnée flottante verticale du pixel (axe Y).
	       z, ///< Coordonnée flottante de profondeur du pixel (axe Z).
	       v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};

/// <summary>Structure représentant des coordonnées (x, y).</summary>
struct Coord {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y; ///< Coordonnée entière verticale du pixel (axe Y).
	[[nodiscard]] size_t GetSize() const { return x * y; }
	// Surcharges arithmétiques avec un autre Coord3D
	Coord operator+(const Coord& other) const { return {.x = x + other.x, .y = y + other.y}; }
	Coord operator-(const Coord& other) const { return {.x = x - other.x, .y = y - other.y}; }
	Coord operator*(const Coord& other) const { return {.x = x * other.x, .y = y * other.y}; }
	Coord operator/(const Coord& other) const { return {.x = x / other.x, .y = y / other.y}; }

	// Surcharges arithmétiques avec un scalaire
	Coord operator+(const size_t scalar) const { return {.x = x + scalar, .y = y + scalar}; }
	Coord operator-(const size_t scalar) const { return {.x = x - scalar, .y = y - scalar}; }
	Coord operator*(const size_t scalar) const { return {.x = x * scalar, .y = y * scalar}; }
	Coord operator/(const size_t scalar) const { return {.x = x / scalar, .y = y / scalar}; }
};

/// <summary>Structure représentant des coordonnées (x, y).</summary>
struct Vec2 {
	double x, ///< Coordonnée flottante horizontale du pixel (axe X).
	       y; ///< Coordonnée flottante verticale du pixel (axe Y).

	// Surcharges arithmétiques avec un autre Vec3
	Vec2 operator+(const Vec2& other) const { return {.x = x + other.x, .y = y + other.y}; }
	Vec2 operator-(const Vec2& other) const { return {.x = x - other.x, .y = y - other.y}; }
	Vec2 operator*(const Vec2& other) const { return {.x = x * other.x, .y = y * other.y}; }
	Vec2 operator/(const Vec2& other) const { return {.x = x / other.x, .y = y / other.y}; }

	// Surcharges arithmétiques avec un scalaire
	Vec2 operator+(const double scalar) const { return {.x = x + scalar, .y = y + scalar}; }
	Vec2 operator-(const double scalar) const { return {.x = x - scalar, .y = y - scalar}; }
	Vec2 operator*(const double scalar) const { return {.x = x * scalar, .y = y * scalar}; }
	Vec2 operator/(const double scalar) const { return {.x = x / scalar, .y = y / scalar}; }
};

/// <summary>Structure représentant des coordonnées (x, y, z).</summary>
struct Coord3D {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y, ///< Coordonnée entière verticale du pixel (axe Y).
	       z; ///< Coordonnée entière de profondeur du pixel (axe Z).

	[[nodiscard]] size_t GetPlaneSize() const { return x * y; }
	[[nodiscard]] size_t GetSize() const { return x * y * z; }
	[[nodiscard]] size_t GetIndex(const Coord3D& size) const { return (z * size.y + y) * size.x + x; }
	[[nodiscard]] Coord3D GetMid() const { return Coord3D{.x = x / 2, .y = y / 2, .z = z / 2}; }

	// Surcharges arithmétiques avec un autre Coord3D
	Coord3D operator+(const Coord3D& other) const { return {.x = x + other.x, .y = y + other.y, .z = z + other.z}; }
	Coord3D operator-(const Coord3D& other) const { return {.x = x - other.x, .y = y - other.y, .z = z - other.z}; }
	Coord3D operator*(const Coord3D& other) const { return {.x = x * other.x, .y = y * other.y, .z = z * other.z}; }
	Coord3D operator/(const Coord3D& other) const { return {.x = x / other.x, .y = y / other.y, .z = z / other.z}; }

	// Surcharges arithmétiques avec un scalaire
	Coord3D operator+(const size_t scalar) const { return {.x = x + scalar, .y = y + scalar, .z = z + scalar}; }
	Coord3D operator-(const size_t scalar) const { return {.x = x - scalar, .y = y - scalar, .z = z - scalar}; }
	Coord3D operator*(const size_t scalar) const { return {.x = x * scalar, .y = y * scalar, .z = z * scalar}; }
	Coord3D operator/(const size_t scalar) const { return {.x = x / scalar, .y = y / scalar, .z = z / scalar}; }
};

/// <summary>Structure représentant des coordonnées (x, y, z).</summary>
struct Vec3 {
	double x, ///< Coordonnée flottante horizontale du pixel (axe X).
	       y, ///< Coordonnée flottante verticale du pixel (axe Y).
	       z; ///< Coordonnée flottante de profondeur du pixel (axe Z).

	// Surcharges arithmétiques avec un autre Vec3
	Vec3 operator+(const Vec3& other) const { return {.x = x + other.x, .y = y + other.y, .z = z + other.z}; }
	Vec3 operator-(const Vec3& other) const { return {.x = x - other.x, .y = y - other.y, .z = z - other.z}; }
	Vec3 operator*(const Vec3& other) const { return {.x = x * other.x, .y = y * other.y, .z = z * other.z}; }
	Vec3 operator/(const Vec3& other) const { return {.x = x / other.x, .y = y / other.y, .z = z / other.z}; }

	// Surcharges arithmétiques avec un scalaire
	Vec3 operator+(const double scalar) const { return {.x = x + scalar, .y = y + scalar, .z = z + scalar}; }
	Vec3 operator-(const double scalar) const { return {.x = x - scalar, .y = y - scalar, .z = z - scalar}; }
	Vec3 operator*(const double scalar) const { return {.x = x * scalar, .y = y * scalar, .z = z * scalar}; }
	Vec3 operator/(const double scalar) const { return {.x = x / scalar, .y = y / scalar, .z = z / scalar}; }
};


//-------------------------------------------------------------------------------------------------
//------------------------------------------ Conversions ------------------------------------------
//-------------------------------------------------------------------------------------------------
/// <summary>Récupère les coordonnées depuis un index et les dimensions du volume.</summary>
/// <param name="index">Index dans le volume.</param>
/// <param name="size">Dimensions du volume.</param>
/// <returns> Un vecteur de position</returns>
inline Coord3D GetCoord3D(size_t index, const Coord3D& size)
{
	Coord3D res{};
	res.x = index % size.x;
	index /= size.x;
	res.y = index % size.y;
	index /= size.y;
	res.z = index;
	return res;
}

/// <summary>Cast des Coordonnées (entières) vers les Coordonnées (flottantes).</summary>
/// <param name="c">Coordonnées entières.</param>
/// <returns> Coordonnées flottantes.</returns>
inline Vec3 CoordToVec(const Coord3D& c) { return Vec3{.x = static_cast<double>(c.x), .y = static_cast<double>(c.y), .z = static_cast<double>(c.z)}; }

/// <summary>Cast des Coordonnées (entières) vers les Coordonnées (flottantes).</summary>
/// <param name="c">Coordonnées entières.</param>
/// <returns> Coordonnées flottantes.</returns>
inline Vec2 CoordToVec(const Coord& c) { return Vec2{.x = static_cast<double>(c.x), .y = static_cast<double>(c.y)}; }

//-------------------------------------------------------------------------------------------------
//------------------------------------------- Surcharges ------------------------------------------
//-------------------------------------------------------------------------------------------------

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Coord3D. Format : (x, y, z) </summary>
inline std::ostream& operator<<(std::ostream& os, const Coord& c)
{
	os << "(" << c.x << ", " << c.y << ")";
	return os;
}

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Coord3D. Format : (x, y, z) </summary>
inline std::ostream& operator<<(std::ostream& os, const Vec2& c)
{
	os << "(" << c.x << ", " << c.y << ")";
	return os;
}

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Coord3D. Format : (x, y, z) </summary>
inline std::ostream& operator<<(std::ostream& os, const Coord3D& c)
{
	os << "(" << c.x << ", " << c.y << ", " << c.z << ")";
	return os;
}

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Coord3D. Format : (x, y, z) </summary>
inline std::ostream& operator<<(std::ostream& os, const Vec3& c)
{
	os << "(" << c.x << ", " << c.y << ", " << c.z << ")";
	return os;
}
} // namespace palm_tracer
