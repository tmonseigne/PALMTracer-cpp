///-------------------------------------------------------------------------------------------------
///
/// \file Interpolation.hpp
/// \brief Interpolation cubique 3D avec support de deux modes :
/// - Uniform scaling (facteur entier par axe)
/// - Interpolation ponctuelle (coordonnées flottantes X, Y, Z)
///
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <vector>
#include "Structs.hpp"


namespace palm_tracer {
/// <summary> Effectue une interpolation trilineaire dans un volume 3D. </summary>
/// <param name="volume">Vecteur linéarisé représentant le volume 3D (X-Y-Z) en mémoire row-major.</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
/// <param name="pos">Position flottante (x, y, z) où calculer l'interpolation.</param>
/// <returns> Valeur interpolée par interpolation linéaire trivariée au point donné.
/// Renvoie 0.0 si la position est en dehors des bornes du volume. </returns>
double LinearInterpolation3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos);

/// <summary>Applique un zoom isotrope sur un volume 3D avec interpolation linéaire.</summary>
/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
/// <param name="level">Facteur de zoom entier (> 0). Un niveau 1 insère un point entre chaque valeur initiale, un niveau 2 en insère 3, etc.
/// Soit `2^factor - 1` points interpolés entre chaque voxel.</param>
/// <param name="newSize">Taille du volume de sortie, calculée comme :
/// - `(size.x - 1) * 2^level + 1`
/// - `(size.y - 1) * 2^level + 1`
/// - `(size.z - 1) * 2^level + 1`
/// Ce comportement imite fidèlement le comportement de MATLAB `interp3(M, level, 'cubic')`.</param>
/// <returns>Volume agrandi sous forme linéarisée avec interpolation linéaire.</returns>
std::vector<double> LinearZoom3D(const std::vector<double>& volume, const Coord3D& size, size_t level, Coord3D& newSize);

/// <summary> Applique un décalage (shift) flottant sur tout le volume avec interpolation linéaire. </summary>
/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
/// <param name="shift">Décalage flottant (x, y, z).</param>
/// <returns>Volume décalé, même taille que l'entrée.</returns>
std::vector<double> LinearShift3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& shift);

/// <summary>Effectue une interpolation cubique 1D sur quatre valeursBasé sur la méthode de Catmull-Rom (forme cardinale).</summary>
/// <param name="v">Tableau de 4 valeurs consécutives.</param>
/// <param name="x">Position relative entre les points (dans [0, 1]).</param>
/// <returns>Valeur interpolée à la position x.</returns>
double CubicInterpolation1D(const std::array<double, 4>& v, double x);

/// <summary>Interpolation cubique trivariée dans un volume 3D (64 voxels pris en compte, 4x4x4).</summary>
/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
/// <param name="pos">Coordonnées flottantes où interpoler.</param>
/// <returns>Valeur interpolée au point donné, ou 0.0 si hors volume.</returns>
double CubicInterpolation3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos);

/// <summary> Applique un décalage (shift) flottant sur tout le volume avec interpolation cubique. </summary>
/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
/// <param name="shift">Décalage flottant (x, y, z).</param>
/// <returns>Volume décalé, même taille que l'entrée.</returns>
std::vector<double> CubicShift3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& shift);
} // namespace palm_tracer
