///-------------------------------------------------------------------------------------------------
///
/// \file Constants.hpp
/// \brief  Fichier de constantes necessaire aux algorithmes.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 18/04/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

namespace palm_tracer {
//=========================================================
//==================== PALM PROCESSING ====================
//=========================================================
constexpr size_t PLANE_ID            = 1;            ///< Indice du plan d'ondelettes à utiliser pour la détection.
constexpr size_t N_PLANES            = PLANE_ID + 2; ///< Nombre de plans d'ondelettes à utiliser pour la détection.
constexpr size_t INTENSITY_MIN       = 0;            ///< Intensité minimale requise pour valider une détection.
constexpr size_t VOLUME_MIN          = 4;            ///< Volume minimal requis pour valider une détection.
constexpr size_t THRESHOLD_ITERATION = 5;            ///< Nombre d'iteration maximum pour l'auto-threshold.

// Paramètres de la localisation :
// Id, Plane, Index, Channel, X, Y, Z, Integrated Intensity, Sigma X, Sigma Y, Theta,
// MSE XY, MSE Z, Intensity 0, Intensity Offset, Intensity, Surface, Circularity
constexpr size_t N_LOCALIZATIONS = 18; ///< Nombre de paramètres à renvoyer lors de la localization.

constexpr double EPSILON = 1e-7;

//==================================================
//==================== TRACKING ====================
//==================================================
// Paramètres du Tracking :
// Track, Plane, Id, X, Y, Z, Integrated Intensity, Surface
constexpr size_t N_TRACKING = 8; ///< Nombre de paramètres à renvoyer lors du tracking.

// Anciennement "Sigma X", "Sigma Y", "Theta", "Y", "X", "Intensity 0", "Intensity Offset", "MSE Gaussian", "Intensity", "Surface", "Z", "Pair Distance", "Id"
// Maintenant Id, X, Y, Z, Integrated Intensity, Surface
constexpr size_t N_IN_PARAMS = 6; ///< Nombre de paramètres reçu lors du tracking.
constexpr size_t N_IN_PARAMS_2 = 4; ///< Nombre de paramètres reçu lors du tracking.

//=================================================
//==================== SPLINES ====================
//=================================================
constexpr size_t SPLINE_PATCH_SIDE  = 4,  ///< Taille d'un côté du patch spline (nombre de points en X ou Y).
                 SPLINE_PATCH_PLANE = 16, ///< Nombre total de points dans un plan 2D du patch spline : SPLINE_PATCH_SIDE ^ 2.
                 SPLINE_PATCH_SIZE  = 64; ///< Nombre total de coefficients dans un patch spline 3D : SPLINE_PATCH_SIDE ^ 3.
}                                         // namespace palm_tracer
