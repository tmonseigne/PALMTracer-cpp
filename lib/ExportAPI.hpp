#pragma once

// Détection de plateforme et définition de la macro EXPORT_API
#if defined(_WIN32) || defined(_WIN64)
#define EXPORT_API __declspec(dllexport)
#define CALL_CONV __stdcall
#else
#define EXPORT_API __attribute__((visibility("default")))
#define CALL_CONV
#endif
