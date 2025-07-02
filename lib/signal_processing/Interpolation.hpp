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
#include <vector>
#include <cmath>
#include <algorithm>


namespace palm_tracer {
// Fonction d'interpolation cubique 1D (Catmull-Rom ou Hermite, ici simple forme cardinal)
double CubicInterpolation1D(double v0, double v1, double v2, double v3, double x);

// Fonction d'interpolation trilineaire cubique (64 valeurs, 4x4x4)
double CubicInterpolation3D(const std::vector<double>& volume, size_t Z, size_t Y, size_t X, double z, double y, double x);

// Fonction de "zoom" uniforme (interpolation cubique avec facteur entier)
std::vector<double> CubicZoom3D(const std::vector<double>& volume, size_t Z, size_t Y, size_t X, int factor);
} // namespace palm_tracer
