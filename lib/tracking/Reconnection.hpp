///-------------------------------------------------------------------------------------------------
/// \file Reconnection.hpp
/// <summary>Reconnexion des segments séparés par scintillement.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <utility>
#include <vector>
#include "TrackingPoint.hpp"

namespace palm_tracer::tracking::reconnection {
/// <summary>Segment cible et coût lexicographique d’une reconnexion admissible.</summary>
struct Candidate {
	size_t target;                  ///< Indice du segment cible dans le conteneur fourni.
	std::pair<size_t, double> cost; ///< Écart de plan, puis distance spatiale.
};

/// <summary>Coût lexicographique d'une reconnexion admissible : frames minimales, puis distance minimale.</summary>
/// <param name="frameGap">Différence de numéro de frame entre les deux extrémités, strictement positive.</param>
/// <param name="distance">Distance spatiale entre les deux extrémités, dans l'unité de <c>Point::pos_</c>.</param>
/// <returns>Paire \f$(\Delta f,d)\f$ a minimiser lexicographiquement ; ce n'est pas une somme pondérée.</returns>
/// <remarks>Les limites de durée et de distance restent des contraintes d'admissibilité séparées du coût.</remarks>
std::pair<size_t, double> GetCost(size_t frameGap, double distance);

/// <summary> Calcule la distance maximale autorisée pour une reconnexion en fonction du <paramref name="mode"/> choisi, de la durée du scintillement et de
/// la vitesse maximale spécifiée.</summary>
/// <param name="mode">Mode de dispersion des points (<see cref="blinking_mode"/>).</param>
/// <param name="duration">Durée du scintillement (\f$\Delta t\f$, en nombre de plans).</param>
/// <param name="maxSpeed">Vitesse maximale supposée (micromètre/plan).</param>
/// <returns>La distance maximale de recherche (en micromètre) pour ce couple (<paramref name="mode"/>, \f$\Delta t\f$).</returns>
double GetMaxDistance(size_t mode, double duration, double maxSpeed);

/// <summary>Génère et trie les reconnexions admissibles d'un segment, indépendamment des conflits gloutons.</summary>
/// <param name="starts">Indices des segments non vides, triés par plan de début puis indice d'origine.</param>
/// <param name="current">Indice du segment dont la fin est à prolonger.</param>
/// <param name="mode">Mode de calcul de la distance maximale admissible.</param>
/// <param name="maxDuration">Délai maximal inclusif, en plans.</param>
/// <param name="maxSpeed">Paramètre spatial transmis à GetMaxDistance.</param>
/// <param name="candidates">Liste vidée puis remplie ; tous les candidats admissibles sont conservés.</param>
/// <returns>Vrai si la liste est triée ; faux si un coût NaN impose de conserver l'ordre temporel initial.</returns>
/// <remarks>Le tri stable utilise exclusivement GetCost. Les filtres taken et indice après la tête
/// sont appliqués lors des connexions séquentielles. Les appels concurrents utilisent des listes de sortie distinctes.</remarks>
/// <param name="tracks">Segments dans leur ordre d’origine ; modifiés uniquement par la reconnexion.</param>
bool GetCandidates(const std::vector<Track>& tracks, const std::vector<size_t>& starts, size_t current, size_t mode,
                   size_t maxDuration, double maxSpeed, std::vector<Candidate>& candidates);

/// <summary>Reconnecte des trajectoires disjointes lorsqu'elles sont séparées par un scintillement (frames manquantes).</summary>
/// <param name="mode">Mode de dispersion utilisé pour calculer la distance maximale de recherche (<see cref="blinking_mode"/>).</param>
/// <param name="maxDuration">Durée maximale (en nombre de plans) entre la fin d'une trajectoire et le début d'une autre pour qu'une reconnexion soit
/// possible.</param>
/// <param name="maxSpeed">Vitesse maximale (en pixels/plan) supposée pour un point entre deux plans consécutifs.</param>
/// <remarks>
/// La reconnexion est effectuée de manière gloutonne :
/// - priorité au plus petit décalage temporel \f$\Delta t\f$,
/// - puis, à  \f$\Delta t\f$ égal, priorité à la distance spatiale minimale,
/// sous réserve que cette distance soit inférieure au seuil calculé par <c>GetMaxDistance</c>.
/// Un index temporel limite la génération parallèle des candidats à la fenêtre admissible.
/// Chaque liste est triée de façon stable par GetCost, puis les connexions sont décidées séquentiellement
/// dans l'ordre des têtes, en vérifiant les segments déjà pris. Les coûts NaN conservent le parcours historique.
/// Dans un même plan, les indices d'origine sont parcourus dans l'ordre croissant pour préserver les ex aequo.
/// Le classement des candidats admissibles est centralise dans <see cref="GetCost"/>.
/// A coût exactement égal, le premier candidat rencontre est conservé ; l'association globale reste gloutonne.
/// </remarks>
/// <param name="tracks">Segments dans leur ordre d’origine ; modifiés uniquement par la reconnexion.</param>
void Reconnect(std::vector<Track>& tracks, size_t mode, size_t maxDuration, double maxSpeed);
} // namespace palm_tracer::tracking::reconnection
