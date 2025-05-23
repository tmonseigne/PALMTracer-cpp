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

// Paramètres de la segmentation
// Plane, X, Y, Z, Integrated Intensity, Sigma X, Sigma Y, Theta,
// MSE Gaussian, MSE Z, Pair Distance, Intensity 0, Intensity Offset, Intensity, Surface, Circularity
constexpr size_t N_SEGMENTS = 16; ///< Nombre de paramètres à renvoyer.

constexpr double EPSILON = 1e-7;

//==================================================
//==================== TRACKING ====================
//==================================================
}
