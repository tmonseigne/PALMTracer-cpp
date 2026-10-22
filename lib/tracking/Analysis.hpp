///-------------------------------------------------------------------------------------------------
/// \file Analysis.hpp
/// <summary>Analyse des trajectoires : MSD, diffusion instantanée et métriques ajustées.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <vector>
#include "ComputesParams.hpp"
#include "TrackingPoint.hpp"

namespace palm_tracer::tracking::analysis {
/// <summary>Calcule le <b>Mean Squared Displacement (MSD)</b>, en considérant tous les décalages temporels (écart de plan) possibles sur une trajectoire.
/// - Soit une trajectoire de longueur \f$N\f$ : positions successives \f$P_1, P_2, …, P_n\f$.
/// - Pour chaque décalage (ou <i>lag</i>) \f$i = 1 … (\mathrm{plan}(P_n) - \mathrm{plan}(P_1))\f$,
///   on calcule la moyenne des distances quadratiques : \f$\mathrm{MSD}(i) = \frac{1}{N_i} \sum_{j=1}^{N_i} \lVert P_{j+i} - P_j \rVert^2\f$
///   avec \f$j\f$ parcourant toutes les paires séparées de \f$i\f$ pas.
/// - Le dernier lag (\f$i = N-1\f$) est donc basé sur un seul terme (\f$\lVert P_{n} - P_1 \rVert^2\f$).
///
/// En cas de blinking, il est possible qu'un lag ne possède aucun ensemble de points le représentant : le MSD est alors fixé à <c>-1</c>.
/// </summary>
/// <param name="track">Trajectoire à analyser.</param>
/// <param name="start">Premier point de la trajectoire à prendre (utile lorsqu'un sous-ensemble de points est sélectionné pour l'ajustement).</param>
/// <param name="end">Dernier point de la trajectoire à prendre (utile lorsqu'un sous-ensemble de points est sélectionné pour l'ajustement).</param>
/// <returns>Un vecteur de longueur \f$L = \mathrm{plan}_{\max}-\mathrm{plan}_{\min}\f$ (ou vide si la trajectoire contient moins de 2 points).</returns>
/// <remarks>
/// - Hypothèse : la trajectoire est triée par plan strictement croissant (aucun doublon de plan).
/// - Unités : <c>Point::GetDistance()</c> fixe l'unité (nm). Appliquer un facteur d'échelle en amont si nécessaire.
/// - Les bacs (lags) sans paires observées sont mis à <c>-1</c>.
/// - Complexité : \f$\mathcal{O}(N^2)\f$ pour \f$N\f$ le nombre de points par trajectoire.
/// </remarks>
[[nodiscard]] std::vector<double> ComputeMSD(const Track& track, size_t start = 0, size_t end = 0);

/// <summary>Calcule en parallèle le MSD de chaque trajectoire et assemble les lignes de résultat.</summary>
/// <param name="tracks">Trajectoires à analyser, conservées dans leur ordre d'origine.</param>
/// <returns>Une ligne par trajectoire : identifiant suivi des MSD par lag, avec -1 pour les lags sans observations.
/// La ligne est vide si la trajectoire contient moins de deux points.</returns>
[[nodiscard]] std::vector<std::vector<double>> ComputeAllMSD(const std::vector<Track>& tracks);

inline constexpr double MinimumInstantDiffusion = 1e-6; ///< Plancher historique de diffusion instantanée.

/// <summary>Estime la <b>diffusion instantanée</b> (Instant D) le long d'une trajectoire, via des <b>fenêtres glissantes</b> de longueur
/// <c>ComputesParams::fitLength</c>.
/// Pour chaque fenêtre, on calcule le vecteur MSD local pour les décalages \f$i = 1..\text{fitLength}\f$, puis on ajuste une droite : \f$\mathrm{MSD}(\tau)
/// \approx 2\,d\,D_{\text{inst}}\,\tau + C\f$, où \f$d\f$ est la dimension (2 ou 3) et \f$\tau = i\,\Delta t\f$.
/// </summary>
/// <param name="track">Trajectoire à analyser.</param>
/// <param name="p">Paramètres de configuration des calculs.</param>
/// <returns>Un vecteur de taille \f$\max(0, N - \text{fitLength})\f$ où \f$N\f$ est la longueur de la trajectoire.</returns>
std::vector<double> ComputeInstantDiffusion(const Track& track, const ComputesParams& p);

/// <summary>Calcule des métriques de diffusion à partir d'une trajectoire et réalise un ou plusieurs fittings du MSD.</summary>
/// <param name="track">Trajectoire à analyser.</param>
/// <param name="p">Paramètres de configuration des calculs.</param>
/// <returns>Tableau contenant le résultat des métriques (taille et contenu dépendant du mode d'ajustement).</returns>
std::vector<double> ComputeFitMetrics(const Track& track, const ComputesParams& p);

/// <summary>Retourne le nombre d'éléments dans le tableau final attendu pour le <see cref="fit_mode"/> courant.</summary>
/// <returns>
/// \f[
/// \begin{cases}
/// 8 & \text{pour } \mathtt{fit\_mode::linear} \quad & ([L, I, D_a, D_b, D_{err}, a, b, err])\\[2pt]
/// 9 & \text{pour } \mathtt{fit\_mode::power} \quad & ([L, I, D_a, D_b, D_{err}, α, a, err, \mu(\text{dist})])\\[2pt]
/// 10 & \text{pour } \mathtt{fit\_mode::exponential} \quad & ([L, I, D_a, D_b, D_{err}, a, \tau, b, err, \sqrt{\frac{3}{4}\times(b+a)}])
/// \end{cases}
/// \f]
/// </returns>
/// <remarks>En cas de <paramref name="mode"/> non reconnu, retourne <c>0</c> (à éviter en pratique).</remarks>
size_t GetFitResultSize(const fit_mode mode);

/// <summary>Calcule en parallèle la diffusion instantanée de chaque trajectoire et assemble les lignes de résultat.</summary>
/// <param name="tracks">Trajectoires à analyser, conservées dans leur ordre d'origine.</param>
/// <param name="p">Paramètres de configuration des calculs.</param>
/// <returns>Une ligne par trajectoire : identifiant suivi des diffusions instantanées.
/// La ligne est vide si aucune fenêtre de diffusion ne peut être calculée.</returns>
[[nodiscard]] std::vector<std::vector<double>> ComputeAllInstantDiffusion(const std::vector<Track>& tracks, const ComputesParams& p);

/// <summary>Calcule en parallèle les métriques ajustées de chaque trajectoire et assemble les lignes de résultat.</summary>
/// <param name="tracks">Trajectoires non vides à analyser, conservées dans leur ordre d'origine.</param>
/// <param name="p">Paramètres de configuration des calculs.</param>
/// <returns>Une ligne par trajectoire : identifiant suivi des métriques du mode choisi.
/// La ligne est vide si l'ajustement ne produit aucune métrique.</returns>
[[nodiscard]] std::vector<std::vector<double>> ComputeAllFitMetrics(const std::vector<Track>& tracks, const ComputesParams& p);
} // namespace palm_tracer::tracking::analysis
