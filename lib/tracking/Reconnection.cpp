///-------------------------------------------------------------------------------------------------
/// \file Reconnection.cpp
/// <summary>Classement des reconnexions et recomposition des trajectoires.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#include "Reconnection.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <ranges>

namespace palm_tracer::tracking::reconnection {
//-------------------------------------------------------------------------------------------------
std::pair<size_t, double> GetCost(const size_t frameGap, const double distance)
{
	// Choix actuel : priorité stricte à la frame la plus proche, puis à la distance à frame égale.
	// La comparaison lexicographique permet à un candidat plus précoce de gagner même s'il est plus éloigné.
	// Variante future possible : C = wt * (frameGap / Tref) + wd * (distance / Dref)^2.
	// Un coût scalaire pourrait être retourné sous la forme {0, C} pour supprimer la priorité stricte au temps.
	// On pourrait ajouter wi * log((I_next + epsilon) / (I_previous + epsilon))^2 pour la continuite d'intensité,
	// en transmettant aussi les intensités. Poids et échelles Tref/Dref seraient à valider physiquement.
	// Ces variantes ne sont pas actives : contraintes de recherche et arbitrage global glouton restent séparés.
	return {frameGap, distance};
}

//-------------------------------------------------------------------------------------------------
double GetMaxDistance(const size_t mode, const double duration, const double maxSpeed)
{
	switch (mode) {
		case 0: return maxSpeed;                                // Immobile
		default: case 1: return maxSpeed * std::sqrt(duration); // Diffusive
		case 2: return maxSpeed * duration;                     // Linear
	}
}

//-------------------------------------------------------------------------------------------------
bool GetCandidates(const std::vector<Track>& tracks, const std::vector<size_t>& starts, const size_t current, const size_t mode,
                   const size_t maxDuration, const double maxSpeed, std::vector<Candidate>& candidates)
{
	candidates.clear();
	if (tracks[current].points.empty()) { return true; }
	const Point& Last = tracks[current].points.back();
	// Premier début strictement postérieur.
	const auto Begin = std::ranges::upper_bound(starts, Last.plane_, {}, [&](const size_t j) { return tracks[j].points.front().plane_; });
	bool ordered     = true;
	for (auto it = Begin; it != starts.end(); ++it) {
		const size_t J     = *it;
		const Point& First = tracks[J].points.front();
		const size_t Gap   = First.plane_ - Last.plane_; // Soustraction sûre ; aucune addition susceptible de déborder.
		if (Gap > maxDuration) { break; }                // L'index est trié : tous les débuts suivants dépassent aussi la fenêtre.
		if (J == current) { continue; }
		const double Distance = First.GetDistance(Last);
		if (Distance > GetMaxDistance(mode, static_cast<double>(Gap), maxSpeed)) { continue; }
		const auto Cost = GetCost(Gap, Distance); // Règle de classement centralisée, indépendante de l'admissibilité.
		candidates.push_back({.target = J, .cost = Cost});
		ordered = ordered && !std::isnan(Cost.second);
	}
	// Les coûts égaux gardent leur ordre temporel puis leur indice d'origine.
	if (ordered) { std::ranges::stable_sort(candidates, {}, &Candidate::cost); }

	return ordered;
}

//-------------------------------------------------------------------------------------------------
void Reconnect(std::vector<Track>& tracks, const size_t mode, const size_t maxDuration, const double maxSpeed)
{
	if (tracks.empty()) { return; }
	const size_t N = tracks.size();
	std::vector taken(N, false);             // True si l'index appartient déjà à une chaîne (hors tête)
	std::vector<std::vector<size_t>> chains; // Liste des chaînes d'indices (ex: { [2,4,5], [3], [6,7], ... })
	chains.reserve(N);                       // Taille maximale possible

	// Index temporel indépendant de l'ordre des têtes : seuls les débuts réellement occupés sont stockés.
	std::vector<size_t> starts; // Indices dans tracks : les trajectoires et l’ordre des têtes restent inchangés.
	starts.reserve(N);
	for (size_t j = 0; j < N; ++j) { if (!tracks[j].points.empty()) { starts.push_back(j); } }
	std::ranges::sort(starts, [&](const size_t a, const size_t b) {
		const size_t A = tracks[a].points.front().plane_, B = tracks[b].points.front().plane_;
		return A != B ? A < B : a < b; // Tri par début ; à plan égal, conserve l’ordre de parcours historique.
	});

	// Chaque itération écrit sa propre liste ; tracks et l'index temporel restent en lecture seule jusqu'à la fin de l'équipe.
	std::vector<std::vector<Candidate>> candidates(N);
	std::vector<unsigned char> ordered(N, 1); // Un octet par liste : vector<bool> partagerait des bits entre threads.
	#pragma omp parallel for schedule(dynamic)
	for (ptrdiff_t current = 0; current < static_cast<ptrdiff_t>(N); ++current) {
		const auto Index = static_cast<size_t>(current);
		ordered[Index]   = static_cast<unsigned char>(GetCandidates(tracks, starts, Index, mode, maxDuration, maxSpeed, candidates[Index]));
	}

	// --- Construction des chaînes d'indices de connexion ---
	for (size_t head = 0; head < N; ++head) {          // Trajectoire à tester
		if (tracks[head].points.empty()) { continue; } // Au cas où, mais théoriquement impossible.
		if (taken[head]) { continue; }                 // Trajectoire déjà absorbée par une chaîne précédente
		std::vector<size_t> chain;                     // Nouvelle Trajectoire
		chain.push_back(head);                         // Ajout de la trajectoire initiale

		// On tente d'allonger la chaîne de façon gloutonne
		size_t current = head; // Segment à prolonger ; head reste la tête de toute la chaîne.
		while (true) {
			size_t bJ     = std::numeric_limits<size_t>::max(); // Indice du meilleur segment encore disponible.
			auto bestCost = GetCost(std::numeric_limits<size_t>::max(), std::numeric_limits<double>::max());
			for (const auto& C : candidates[current]) {
				if (C.target <= head || taken[C.target]) { continue; } // Les conflits restent résolus dans l'ordre des têtes.
				if (C.cost >= bestCost) { continue; }
				bJ       = C.target;
				bestCost = C.cost;
				if (ordered[current] != 0) { break; } // Le premier candidat libre est optimal dans une liste triée.
				// Un coût NaN impose le parcours complet historique, sans lui attribuer artificiellement une priorité.
			}

			if (bJ == std::numeric_limits<size_t>::max()) { break; } // Aucune trajectoire de trouvé, on sort du while
			chain.push_back(bJ);                                     // Ajout de la trajectoire à la liste
			taken[bJ] = true;                                        // Réserve ce segment pour cette chaîne et empêche sa réutilisation.
			current   = bJ;                                          // Mise à jour de la trajectoire courante (qui est la dernière de la liste)
		}
		chains.emplace_back(std::move(chain)); // Ajout de la chaîne de trajectoire au tableau final.
	}

	// --- Recomposition finale dans un nouveau conteneur, par move-concat ---
	std::vector<Track> new_tracks;     // Création d'une liste de trajectoires temporaire
	new_tracks.reserve(chains.size()); // On a notre nombre de trajectoires finales
	size_t track_id = 1;
	for (const auto& Ch : chains) {
		// Concatène tous les points des indices de la chaîne dans un seul track
		Track acc;
		acc.id = track_id++;
		// acc.id = tracks[Ch.front()].id; si l'on veut garder les numéro précédents, mais il y aura des trou dnas les numéros de trajectoires.

		// (option : size hint) somme des tailles pour limiter les reallocs
		size_t totalPts = 0;
		for (const size_t Idx : Ch) { totalPts += tracks[Idx].points.size(); }
		acc.points.reserve(totalPts);

		for (const size_t Idx : Ch) {
			auto& tr = tracks[Idx];
			acc.points.insert(acc.points.end(), std::make_move_iterator(tr.points.begin()), std::make_move_iterator(tr.points.end()));
			tr.points.clear(); // libère les références, pas la capacité (sans importance ici)
		}
		new_tracks.emplace_back(std::move(acc));
	}
	tracks.swap(new_tracks); // tracks devient le résultat reconstruit (complexité : O(1))
}
} // namespace palm_tracer::tracking::reconnection
