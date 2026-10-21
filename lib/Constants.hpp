///-------------------------------------------------------------------------------------------------
///
/// \file Constants.hpp
/// <summary>Fichier de constantes nécessaire aux algorithmes.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 18/04/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <stddef.h>

namespace palm_tracer {
//=========================================================
//==================== PALM PROCESSING ====================
//=========================================================
constexpr size_t PLANE_ID            = 1;            ///< Indice du plan d'ondelettes à utiliser pour la détection.
constexpr size_t N_PLANES            = PLANE_ID + 2; ///< Nombre de plans d'ondelettes à utiliser pour la détection.
constexpr size_t INTENSITY_MIN       = 0;            ///< Intensité minimale requise pour valider une détection.
constexpr size_t VOLUME_MIN          = 4;            ///< Volume minimal requis pour valider une détection.
constexpr size_t THRESHOLD_ITERATION = 5;            ///< Nombre d'itération maximum pour l'auto-threshold.

// Paramètres de la localisation :
// Id, Plane, Index, Channel, X, Y, Z, Integrated Intensity, Sigma X, Sigma Y, Theta,
// MSE XY, MSE Z, Intensity 0, Intensity Offset, Intensity, Surface, Circularity

/// <summary>Indices des colonnes d'une localisation dans le tableau DLL.</summary>
namespace localization_columns {
constexpr size_t Id           = 0;  ///< Identifiant global, attribué à partir de 1.
constexpr size_t Plane        = 1;  ///< Numéro du plan, à partir de 1.
constexpr size_t Index        = 2;  ///< Identifiant local au plan, à partir de 1.
constexpr size_t Channel      = 3;  ///< Canal ; -1 en l'absence d'attribution.
constexpr size_t X            = 4;  ///< Position X en pixels.
constexpr size_t Y            = 5;  ///< Position Y en pixels.
constexpr size_t Z            = 6;  ///< Position axiale dans l'unité du modèle.
constexpr size_t Intensity    = 7;  ///< Intensité intégrée selon le modèle.
constexpr size_t SigmaX       = 8;  ///< Écart type selon X en pixels.
constexpr size_t SigmaY       = 9;  ///< Écart type selon Y en pixels.
constexpr size_t Theta        = 10; ///< Angle du modèle gaussien en radians.
constexpr size_t QualityXY    = 11; ///< Indicateur XY ; colonne historique MSE XY.
constexpr size_t QualityZ     = 12; ///< Indicateur axial ; colonne historique MSE Z.
constexpr size_t Amplitude    = 13; ///< Amplitude selon le modèle.
constexpr size_t Offset       = 14; ///< Fond selon le modèle.
constexpr size_t SegIntensity = 15; ///< Intensité de la segmentation.
constexpr size_t Surface      = 16; ///< Surface de la segmentation en pixels carrés.
constexpr size_t Circularity  = 17; ///< Rapport du plus petit au plus grand écart type, ou 1 si indéfini.
constexpr size_t Count        = 18; ///< Nombre de colonnes d'une localisation.
}

constexpr double EPSILON = 1e-7;

//==================================================
//==================== TRACKING ====================
//==================================================
// Paramètres du Tracking :
// Track, Plane, Id, X, Y, Z, Integrated Intensity, Surface

/// <summary>Indices des colonnes des trajectoires échangées avec la DLL.</summary>
namespace tracking_columns {
constexpr size_t Track     = 0; ///< Identifiant de trajectoire.
constexpr size_t Plane     = 1; ///< Numéro du plan.
constexpr size_t Id        = 2; ///< Identifiant de localisation.
constexpr size_t X         = 3; ///< Position X en pixels.
constexpr size_t Y         = 4; ///< Position Y en pixels.
constexpr size_t Z         = 5; ///< Position axiale, conservée sans conversion par le suivi.
constexpr size_t Intensity = 6; ///< Intensité intégrée.
constexpr size_t Surface   = 7; ///< Surface en pixels carrés.
constexpr size_t Count     = 8; ///< Nombre de colonnes d'une trajectoire.
}

/// <summary>Indices des colonnes de l'entrée compacte servant à associer les points.</summary>
namespace tracking_input_columns {
constexpr size_t Id    = 0; ///< Identifiant de localisation.
constexpr size_t Plane = 1; ///< Numéro du plan.
constexpr size_t X     = 2; ///< Position X en pixels.
constexpr size_t Y     = 3; ///< Position Y en pixels.
constexpr size_t Count = 4; ///< Nombre de colonnes de l'entrée compacte.
}

/// <summary>Indices des colonnes de la sortie d'association des points.</summary>
namespace tracking_output_columns {
constexpr size_t Id    = 0; ///< Identifiant de localisation.
constexpr size_t Track = 1; ///< Identifiant de trajectoire.
constexpr size_t Count = 2; ///< Nombre de colonnes de la sortie d'association.
}

//=================================================
//==================== SPLINES ====================
//=================================================
constexpr size_t SPLINE_PATCH_SIDE  = 4,  ///< Taille d'un côté du patch spline (nombre de points en X ou Y).
                 SPLINE_PATCH_PLANE = 16, ///< Nombre total de points dans un plan 2D du patch spline : SPLINE_PATCH_SIDE ^ 2.
                 SPLINE_PATCH_SIZE  = 64; ///< Nombre total de coefficients dans un patch spline 3D : SPLINE_PATCH_SIDE ^ 3.
}                                         // namespace palm_tracer
