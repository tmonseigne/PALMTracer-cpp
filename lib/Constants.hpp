///-------------------------------------------------------------------------------------------------
///
/// \file Constants.hpp
/// \brief  Fichier de constantes necessaire aux algorithmes.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

namespace palm_tracer {
//=========================================================
//==================== PALM PROCESSING ====================
//=========================================================
constexpr size_t PLANE_ID      = 1;            ///< Indice du plan d'ondelettes à utiliser pour la détection.
constexpr size_t N_PLANES      = PLANE_ID + 2; ///< Nombre de plans d'ondelettes à utiliser pour la détection.
constexpr size_t INTENSITY_MIN = 0;            ///< Intensité minimale requise pour valider une détection.
constexpr size_t VOLUME_MIN    = 4;            ///< Volume minimal requis pour valider une détection.

// Paramètres de la localisation :
// Id, Plane, Index, Channel, X, Y, Z, Integrated Intensity,
// Sigma X, Sigma Y, Theta, MSE Gaussian, MSE Z, Pair Distance,
// Intensity 0, Intensity Offset, Intensity, Surface, Circularity
constexpr size_t N_LOCALIZATIONS = 19; ///< Nombre de paramètres à renvoyer lors de la localization.

constexpr double EPSILON = 1e-7;

//==================================================
//==================== TRACKING ====================
//==================================================
// Paramètres du Tracking :
// Track, Plane, Id, X, Y, Z, Integrated Intensity, Surface, Color
constexpr size_t N_TRACKING = 9; ///< Nombre de paramètres à renvoyer lors du tracking.

// Anciennement "Sigma X", "Sigma Y", "Theta", "Y", "X", "Intensity 0", "Intensity Offset", "MSE Gaussian", "Intensity", "Surface", "Z", "Pair Distance", "Id"
// Maintenant Id, X, Y, Z, Integrated Intensity, Pair Distance (color), surface
constexpr size_t N_IN_PARAMS = 7; ///< Nombre de paramètres reçu lors du tracking.
}                                 // namespace palm_tracer                                 // namespace palm_tracer
