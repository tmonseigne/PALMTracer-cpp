///-------------------------------------------------------------------------------------------------
///
/// \file PSF.hpp
/// \brief Fonctions utiles pour le calcul des PSF.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

namespace palm_tracer {
/// <summary>PSF (Point Spread Function) gaussienne isotrope (2D) centrée en (0,0), à sigma fixe.</summary>
/// <param name="offset">Valeur de fond (bruit ou signal de fond).</param>
/// <param name="i0">Intensité maximale de la gaussienne.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigma">Écart-type commun (isotrope).</param>
/// <returns>Valeur de la PSF en (x, y).</returns>
double PSF(double offset, double i0, double x, double y, double sigma);

/// <summary>PSF (Point Spread Function) gaussienne anisotrope 2D avec rotation (theta).</summary>
/// <param name="offset">Valeur de fond (bruit ou signal de fond).</param>
/// <param name="i0">Intensité maximale de la gaussienne.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigmaX">Écart-type selon l'axe X.</param>
/// <param name="sigmaY">Écart-type selon l'axe Y.</param>
/// <param name="theta">Angle de rotation (en radians).</param>
/// <returns>Valeur de la PSF en (x, y).</returns>
double PSF(double offset, double i0, double x, double y, double sigmaX, double sigmaY, double theta);
} // namespace palm_tracer
