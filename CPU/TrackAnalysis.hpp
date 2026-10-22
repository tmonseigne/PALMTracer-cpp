///-------------------------------------------------------------------------------------------------
/// \file TrackAnalysis.hpp
/// <summary>Adaptation des buffers DLL pour l’analyse des trajectoires.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#pragma once

#include <cstdint>

namespace palm_tracer::cpu {
/// <summary>Charge les trajectoires, sélectionne les analyses et sérialise successivement leurs résultats.</summary>
/// <param name="input">Tableau contigu de trajectoires au format historique à huit colonnes.</param>
/// <param name="outputMSD">Sortie irrégulière du MSD ; ignorée si le calcul est désactivé.</param>
/// <param name="outputInd">Sortie irrégulière de diffusion instantanée ; ignorée si le calcul est désactivé.</param>
/// <param name="outputFit">Sortie irrégulière des métriques ; ignorée si le mode vaut zéro.</param>
/// <param name="nRow">Nombre de localisations dans l’entrée.</param>
/// <param name="isMSD">Active le calcul des déplacements quadratiques moyens.</param>
/// <param name="isInd">Active le calcul de diffusion instantanée.</param>
/// <param name="is3D">Active les conventions de diffusion tridimensionnelle.</param>
/// <param name="pixelSize">Échelle des coordonnées XY en unités physiques par pixel.</param>
/// <param name="exposureTime">Intervalle entre deux plans, en secondes.</param>
/// <param name="fitMode">Zéro désactive les fits ; 2 sélectionne la puissance, 3 l’exponentielle, sinon le modèle linéaire.</param>
/// <param name="fitParams">Paramètres empruntés ; le premier élément contient la longueur d’ajustement.</param>
/// <returns>Vrai lorsque l’exécution se termine normalement.</returns>
/// <remarks>Les tampons activés doivent être suffisamment grands pour leurs résultats et la sentinelle finale.
/// Les calculs sont exécutés dans l’ordre historique MSD, diffusion instantanée, puis fits ; chaque tableau temporaire est libéré avant le suivant.
/// L’initialisation globale d’OpenMP reste à la charge de l’export DLL.</remarks>
bool AnalyzeTracks(const double* input, double* outputMSD, double* outputInd, double* outputFit, uint64_t nRow,
                   bool isMSD, bool isInd, bool is3D,
                   double pixelSize, double exposureTime, uint64_t fitMode, const double* fitParams);
} // namespace palm_tracer::cpu
