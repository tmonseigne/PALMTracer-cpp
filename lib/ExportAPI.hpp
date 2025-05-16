///-------------------------------------------------------------------------------------------------
/// 
/// \file ExportAPI.hpp
/// \brief Définition de macros pour l’exportation des symboles de bibliothèque dynamique.
/// 
/// Ce fichier configure la macro `EXPORT_API` en fonction du système de compilation.
/// Elle permet de marquer les fonctions ou classes exportées depuis une DLL (sous Windows) ou
/// d’assurer la visibilité des symboles sur d’autres plateformes (Unix, macOS).
/// Une convention d'appel `CALL_CONV` est également définie pour garantir la compatibilité
/// avec les appels inter-langages ou inter-modules.
/// 
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// 
///-------------------------------------------------------------------------------------------------

#pragma once

// Détection de plateforme et définition de la macro EXPORT_API
#if defined(_WIN32) || defined(_WIN64)
#define EXPORT_API __declspec(dllexport)
#define CALL_CONV __stdcall
#else
#define EXPORT_API __attribute__((visibility("default")))
#define CALL_CONV
#endif
