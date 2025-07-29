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
#include <algorithm>
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
struct Vec2I {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y; ///< Coordonnée entière verticale du pixel (axe Y).

	// Surcharges arithmétiques avec un autre Vec3I
	Vec2I operator+(const Vec2I& other) const { return {.x = x + other.x, .y = y + other.y}; }
	Vec2I operator-(const Vec2I& other) const { return {.x = x - other.x, .y = y - other.y}; }
	Vec2I operator*(const Vec2I& other) const { return {.x = x * other.x, .y = y * other.y}; }
	Vec2I operator/(const Vec2I& other) const { return {.x = x / other.x, .y = y / other.y}; }

	// Surcharges arithmétiques avec un scalaire
	Vec2I operator+(const size_t scalar) const { return {.x = x + scalar, .y = y + scalar}; }
	Vec2I operator-(const size_t scalar) const { return {.x = x - scalar, .y = y - scalar}; }
	Vec2I operator*(const size_t scalar) const { return {.x = x * scalar, .y = y * scalar}; }
	Vec2I operator/(const size_t scalar) const { return {.x = x / scalar, .y = y / scalar}; }

	[[nodiscard]] size_t GetSize() const { return x * y; }
	[[nodiscard]] Vec2I Clamp(const Vec2I& other) const { return {.x = std::clamp<size_t>(x, 0, other.x), .y = std::clamp<size_t>(y, 0, other.y)}; }
};

/// <summary>Structure représentant des coordonnées (x, y).</summary>
struct Vec2D {
	double x, ///< Coordonnée flottante horizontale du pixel (axe X).
	       y; ///< Coordonnée flottante verticale du pixel (axe Y).

	// Surcharges arithmétiques avec un autre Vec3D
	Vec2D operator+(const Vec2D& other) const { return {.x = x + other.x, .y = y + other.y}; }
	Vec2D operator-(const Vec2D& other) const { return {.x = x - other.x, .y = y - other.y}; }
	Vec2D operator*(const Vec2D& other) const { return {.x = x * other.x, .y = y * other.y}; }
	Vec2D operator/(const Vec2D& other) const { return {.x = x / other.x, .y = y / other.y}; }

	// Surcharges arithmétiques avec un scalaire
	Vec2D operator+(const double scalar) const { return {.x = x + scalar, .y = y + scalar}; }
	Vec2D operator-(const double scalar) const { return {.x = x - scalar, .y = y - scalar}; }
	Vec2D operator*(const double scalar) const { return {.x = x * scalar, .y = y * scalar}; }
	Vec2D operator/(const double scalar) const { return {.x = x / scalar, .y = y / scalar}; }
};

/// <summary>Structure représentant des coordonnées (x, y, z).</summary>
struct Vec3I {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y, ///< Coordonnée entière verticale du pixel (axe Y).
	       z; ///< Coordonnée entière de profondeur du pixel (axe Z).

	[[nodiscard]] size_t GetPlaneSize() const { return x * y; }
	[[nodiscard]] size_t GetSize() const { return x * y * z; }
	[[nodiscard]] size_t GetIndex(const Vec3I& size) const { return (z * size.y + y) * size.x + x; }
	[[nodiscard]] Vec3I GetMid() const { return Vec3I{.x = x / 2, .y = y / 2, .z = z / 2}; }

	[[nodiscard]] Vec3I Clamp(const Vec3I& other) const
	{
		return {.x = std::clamp<size_t>(x, 0, other.x), .y = std::clamp<size_t>(y, 0, other.y), .z = std::clamp<size_t>(z, 0, other.z)};
	}

	// Surcharges arithmétiques avec un autre Vec3I
	Vec3I operator+(const Vec3I& other) const { return {.x = x + other.x, .y = y + other.y, .z = z + other.z}; }
	Vec3I operator-(const Vec3I& other) const { return {.x = x - other.x, .y = y - other.y, .z = z - other.z}; }
	Vec3I operator*(const Vec3I& other) const { return {.x = x * other.x, .y = y * other.y, .z = z * other.z}; }
	Vec3I operator/(const Vec3I& other) const { return {.x = x / other.x, .y = y / other.y, .z = z / other.z}; }

	// Surcharges arithmétiques avec un scalaire
	Vec3I operator+(const size_t scalar) const { return {.x = x + scalar, .y = y + scalar, .z = z + scalar}; }
	Vec3I operator-(const size_t scalar) const { return {.x = x - scalar, .y = y - scalar, .z = z - scalar}; }
	Vec3I operator*(const size_t scalar) const { return {.x = x * scalar, .y = y * scalar, .z = z * scalar}; }
	Vec3I operator/(const size_t scalar) const { return {.x = x / scalar, .y = y / scalar, .z = z / scalar}; }

	/// Accès en lecture à une coordonnée par indice (0 → x, 1 → y, 2 → z).
	size_t operator[](const size_t index) const
	{
		switch (index) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: throw std::out_of_range("Vec3D index must be 0, 1 or 2.");
		}
	}

	/// Accès en écriture à une coordonnée par indice (0 → x, 1 → y, 2 → z).
	size_t& operator[](const size_t index)
	{
		switch (index) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: throw std::out_of_range("Vec3D index must be 0, 1 or 2.");
		}
	}
};

/// <summary>Structure représentant des coordonnées (x, y, z).</summary>
struct Vec3D {
	double x, ///< Coordonnée flottante horizontale du pixel (axe X).
	       y, ///< Coordonnée flottante verticale du pixel (axe Y).
	       z; ///< Coordonnée flottante de profondeur du pixel (axe Z).

	// Surcharges arithmétiques avec un autre Vec3D
	Vec3D operator+(const Vec3D& other) const { return {.x = x + other.x, .y = y + other.y, .z = z + other.z}; }
	Vec3D operator-(const Vec3D& other) const { return {.x = x - other.x, .y = y - other.y, .z = z - other.z}; }
	Vec3D operator*(const Vec3D& other) const { return {.x = x * other.x, .y = y * other.y, .z = z * other.z}; }
	Vec3D operator/(const Vec3D& other) const { return {.x = x / other.x, .y = y / other.y, .z = z / other.z}; }

	// Surcharges arithmétiques avec un scalaire
	Vec3D operator+(const double scalar) const { return {.x = x + scalar, .y = y + scalar, .z = z + scalar}; }
	Vec3D operator-(const double scalar) const { return {.x = x - scalar, .y = y - scalar, .z = z - scalar}; }
	Vec3D operator*(const double scalar) const { return {.x = x * scalar, .y = y * scalar, .z = z * scalar}; }
	Vec3D operator/(const double scalar) const { return {.x = x / scalar, .y = y / scalar, .z = z / scalar}; }

	/// Accès en lecture à une coordonnée par indice (0 → x, 1 → y, 2 → z).
	double operator[](const size_t index) const
	{
		switch (index) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: throw std::out_of_range("Vec3D index must be 0, 1 or 2.");
		}
	}

	/// Accès en écriture à une coordonnée par indice (0 → x, 1 → y, 2 → z).
	double& operator[](const size_t index)
	{
		switch (index) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: throw std::out_of_range("Vec3D index must be 0, 1 or 2.");
		}
	}
};


//-------------------------------------------------------------------------------------------------
//------------------------------------------ Conversions ------------------------------------------
//-------------------------------------------------------------------------------------------------
/// <summary>Récupère les coordonnées depuis un index et les dimensions du volume.</summary>
/// <param name="index">Index dans le volume.</param>
/// <param name="size">Dimensions du volume.</param>
/// <returns> Un vecteur de position</returns>
inline Vec3I GetCoord3D(size_t index, const Vec3I& size)
{
	Vec3I res{.x = index % size.x, .y = 0, .z = 0};
	index /= size.x;
	res.y = index % size.y;
	index /= size.y;
	res.z = index;
	return res;
}

/// <summary>Cast des Coordonnées (entières) vers les Coordonnées (flottantes).</summary>
/// <param name="c">Coordonnées entières.</param>
/// <returns> Coordonnées flottantes.</returns>
inline Vec3D CoordToVec(const Vec3I& c) { return Vec3D{.x = static_cast<double>(c.x), .y = static_cast<double>(c.y), .z = static_cast<double>(c.z)}; }

/// <summary>Cast des Coordonnées (entières) vers les Coordonnées (flottantes).</summary>
/// <param name="c">Coordonnées entières.</param>
/// <returns> Coordonnées flottantes.</returns>
inline Vec2D CoordToVec(const Vec2I& c) { return Vec2D{.x = static_cast<double>(c.x), .y = static_cast<double>(c.y)}; }

//-------------------------------------------------------------------------------------------------
//------------------------------------------- Surcharges ------------------------------------------
//-------------------------------------------------------------------------------------------------

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Vec3I. Format : (x, y) </summary>
inline std::ostream& operator<<(std::ostream& os, const Vec2I& c)
{
	os << "(" << c.x << ", " << c.y << ")";
	return os;
}

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Vec3I. Format : (x, y) </summary>
inline std::ostream& operator<<(std::ostream& os, const Vec2D& c)
{
	os << "(" << c.x << ", " << c.y << ")";
	return os;
}

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Vec3I. Format : (x, y, z) </summary>
inline std::ostream& operator<<(std::ostream& os, const Vec3I& c)
{
	os << "(" << c.x << ", " << c.y << ", " << c.z << ")";
	return os;
}

/// <summary>Surcharge de l'opérateur d'affichage pour la structure Vec3I. Format : (x, y, z) </summary>
inline std::ostream& operator<<(std::ostream& os, const Vec3D& c)
{
	os << "(" << c.x << ", " << c.y << ", " << c.z << ")";
	return os;
}
} // namespace palm_tracer
