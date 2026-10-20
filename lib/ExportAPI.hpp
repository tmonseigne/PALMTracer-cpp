///-------------------------------------------------------------------------------------------------
///
/// \file ExportAPI.hpp
/// <summary>Définition de macros pour l'exportation des symboles de bibliothèque dynamique.
///
/// <c>PT_EXPORT</c> assure la visibilité des symboles exportés selon la plateforme.
/// <c>PT_EXTERN_C</c> impose la liaison C, et <c>PT_API</c> combine ces deux marqueurs pour les points d’entrée de la DLL.
/// <c>PT_CALL</c> est actuellement vide et laisse la convention d’appel par défaut du compilateur.
///
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

// Détection de plateforme et définition de la visibilité des symboles.
#if defined(_WIN32) || defined(_WIN64)
#define PT_EXPORT __declspec(dllexport)
#else
#define PT_EXPORT __attribute__((visibility("default")))
#endif

// Marqueurs de convention d’appel et de liaison pour les points d’entrée publics.
#define PT_CALL
#define PT_EXTERN_C extern "C"
#define PT_API PT_EXTERN_C PT_EXPORT
