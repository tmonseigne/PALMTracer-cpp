///-------------------------------------------------------------------------------------------------
///
/// \file Tracking.hpp
/// \brief Définition de la classe <see cref="Tracking"/> qui construit des trajectoires à partir d'une liste de points détectés image par image.
/// Elle lit une liste de points (par frame), établit des connexions temporelles (voisinage avant/après), optimise des liaisons par itérations pseudo-aléatoires, puis exporte le résultat sous forme tabulaire.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 18/04/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once


#include "Parameters.hpp"
#include "TimeNodes.hpp"


namespace palm_tracer::tracking {
/// <summary>Classe réalisant le suivi (tracking) de structures ponctuelles au fil du temps.
/// Les points sont regroupés par instant (<c>TimeNodes</c>). Le suivi construit des liaisons <c>Node::past_</c>, <c>Node::future_</c> entre <c>Node</c>.
/// Le résultat final est un tableau 2D contigu (Track, Plan, Id, X, Y, Z, Intensité, Surface).</summary>
class Tracking {
public:
	/// <summary>Constructeur de la classe <see cref="Tracking"/>.</summary>
	/// <param name="points">Tableau d'entrée contigu (par instants) décrivant les détections, au format <see cref="N_IN_PARAMS"/> par point.
	/// Chaque frame se termine par une ligne sentinelle (Surface = 0 ou -1) indiquant la fin de la frame.</param>
	/// <param name="tracks">Tableau de sortie contigu (pré-alloué) pour les trajectoires (8 colonnes par point).</param>
	/// <param name="minLife">Durée minimale (en nombre de frames) pour qu'une structure soit considérée comme une trajectoire valide.</param>
	/// <param name="decrease">Facteur de décroissance appliqué à l'erreur au dernier pas d'optimisation.</param>
	/// <param name="ntimes">Nombre total d'instants (frames) dans l'entrée.</param>
	/// <param name="params">Paramètres partagés (statistiques et constantes) mis à jour au cours du suivi.</param>
	explicit Tracking(const double* points, double* tracks, const size_t minLife, const double decrease, const size_t ntimes, Parameters* params)
		: points_(points), tracks_(tracks), decrease_(decrease), minLife_(minLife), nbTimes_(ntimes), params_(params) { }


	/// <summary>Destructeur par défaut de la classe <see cref="Tracking"/> : libère le conteneur de frames.</summary>
	~Tracking() { timeNodes_.clear(); }


	/// <summary>Lance le processus complet de tracking.</summary>
	/// <returns>Nombre total de valeurs écrites dans <see cref="tracks_"/> (soit 8 × nombre de points exportés).</returns>
	/// <remarks>
	/// Étapes :
	/// -# Lecture/partition de la liste de points par frame (<c>readPointsList()</c>).
	/// -# Création des voisinages temporels (<c>linkNeighborsBeforeAndAfter()</c>).
	/// -# Tri par intensité des voisinages (<c>sortByIntensityAllNeighbors()</c>).
	/// -# Boucle d'optimisation en <c>Parameters::N_STEPS</c> pas, avec <c>Parameters::N_ITERATIONS</c> itérations par pas :
	///		- initialisation des connexions, estimation des constantes (<c>computeConstants()</c>),
	///		- mises à jour successives des constantes (<c>updateConstants()</c>),
	///		- tentatives aléatoires de modification (<c>TimeNodes::TryRandomChange()</c>).
	/// -# Étiquetage des composantes suivies (<c>setIndexToTrackedStructures()</c>, <c>setColorToTrackedStructures()</c>).
	/// -# Export ASCII/contigu (<c>printConnectionsASCII_trc()</c>).
	/// </remarks>
	size_t Processing();


	size_t nTraj_ = 0; ///< Nombre de trajectoires détectées/étiquetées (mis à jour lors de l'export).

protected:
	// --- Données d'entrée/sortie ---
	const double* points_; ///< Pointeur vers le tableau d'entrée (détections).
	double* tracks_;       ///< Pointeur vers le tableau de sortie (trajectoires).


	// --- Paramètres de contrôle ---
	double decrease_;    ///< Facteur de décroissance de l'erreur au dernier pas.
	size_t minLife_;     ///< Durée minimale d'une trajectoire valide (en frames).
	size_t nbTimes_;     ///< Nombre total de frames.
	Parameters* params_; ///< Paramètres/statistiques partagés.


	// --- Représentation temporelle ---
	std::vector<TimeNodes> timeNodes_; ///< Liste des frames et de leurs nœuds (détections).


	// --- Sous-routines internes ---
	/// <summary>Met à jour les constantes statistiques utilisées dans les coûts (écarts-types et covariance) à partir des liaisons actuelles entre nœuds consécutifs.</summary>
	/// <remarks>
	/// Parcourt toutes les paires (<c>node</c>, <c>node.future_</c>) existantes pour agréger :
	/// - \f$\sigma_{\text{dist}} = \operatorname{std}(\text{Distance}(t\!\to\!t{+}1))\f$
	/// - \f$\sigma_{I} = \operatorname{std}(I_{t} - I_{t{+}1})\f$
	/// - \f$\sigma_{\mu} = \operatorname{std}(\mu_{t} - \mu_{t{+}1})\quad (\mu = I/\text{surf})\f$
	/// - \f$\operatorname{cov}(I,\mu)\f$
	///
	/// Met ensuite à jour ces valeurs dans <see cref="params_"/> et propage vers chaque nœud via <c>Node::UpdateConstants()</c>.\n
	/// Effets de bord : modifie <see cref="params_"/> et le champ « constantes » des nœuds.\n
	/// Complexité : \f$\mathcal{O}(N)\f$ avec \f$N\f$ le nombre total de liaisons.
	/// </remarks>
	void updateConstants();


	/// <summary>Estime par convergence les constantes statistiques initiales
	/// (écarts-types distance/intensité/mean) avant optimisation des liaisons.</summary>
	/// <remarks>
	/// Exécute plusieurs passes d'agrégation entre frames successives (sans dépendre d'un appariement final) pour obtenir des valeurs stables,
	/// puis applique une pondération empirique sur \f$\sigma_{\text{dist}}\f$ afin d'augmenter l'influence de la distance dans le coût.\n
	/// Effets de bord : initialise <c>params_->stdIntensity_</c>, <c>params_->stdMean_</c>, <c>params_->stdDistance_</c>, puis appelle <c>Node::UpdateConstants()</c> pour chaque nœud.\n
	/// Complexité : \f$\mathcal{O}(N)\f$ sur l'ensemble des nœuds.
	/// </remarks>
	void computeConstants();

	/// <summary>Crée les voisinages temporels entre frames consécutives (<c>t → t + 1</c>), permettant d'évaluer les candidats de connexion.</summary>
	/// <remarks>
	/// Pour chaque frame <c>t</c>, appelle <c><see cref="TimeNodes::LinkNeighborsBeforeAndAfter"/>(next)</c> avec la frame <c>t + 1</c>.\n
	/// Les critères de voisinage (fenêtre spatiale, contraintes) sont définis par <c>TimeNodes</c>/<c>Parameters</c>.\n
	/// Effets de bord : remplit les listes de voisins « before/after » dans chaque nœud.\n
	/// Complexité : dépend de l'implémentation (typiquement \f$\mathcal{O}(n_t \times n_{t+1})\f$ sans structure d'indexation).
	/// </remarks>
	void linkNeighborsBeforeAndAfter() { for (size_t t = 0; t < nbTimes_ - 1; ++t) { timeNodes_[t].LinkNeighborsBeforeAndAfter(timeNodes_[t + 1]); } }

	/// <summary>Trie, pour chaque frame, les voisins de chaque nœud par intensité décroissante.</summary>
	/// <remarks>
	/// Ce tri stabilise certains choix gloutons et accélère les heuristiques (accès au meilleur voisin d'abord).\n
	/// Effets de bord : réordonne in-place les listes de voisins dans <see cref="timeNodes_"/>.\n
	/// Complexité : \f$\sum_t \mathcal{O}\left(\sum_{i\in t} \times d_i \times \log(d_i)\right)\f$ où \f$d_i\f$ est le degré de voisinage du nœud \f$i\f$.
	/// </remarks>
	void sortByIntensityAllNeighbors() { for (size_t t = 0; t < nbTimes_; ++t) { timeNodes_[t].SortByIntensityAllNeighbors(); } }

	/// <summary>Lit la liste contiguë des points d'entrée et construit la représentation temporelle <see cref="timeNodes_"/>.</summary>
	/// <remarks>
	/// Hypothèses d'entrée : blocs de points consécutifs par frame, terminés par une ligne sentinelle (Surface = 0 ou −1).\n
	/// Pour chaque point : initialise <c>Node</c> (id, x, y, z, intensity, surf, mean = intensity/surf), met à jour <c>params_->maxX_/maxY_/maxZ_</c>,
	/// puis trie les nœuds de la frame par intensité décroissante.\n
	/// Effets de bord : remplit <see cref="timeNodes_"/>, met à jour <see cref="params_"/>, pré-calcule la distance aux bords via <c>Node::ComputeBorderDistance()</c>.\n
	/// Complexité : \f$\mathcal{O}(N)\f$ pour la lecture + coût du tri par frame.
	/// </remarks>
	void readPointsList();

	/// <summary>Exporte les chaînes suivies sous forme tabulaire contiguë (ASCII-like) dans <see cref="tracks_"/>.</summary>
	/// <returns>Nombre total de valeurs écrites (8 × nombre de points exportés).</returns>
	/// <remarks>
	/// Pour chaque tête de trajectoire (<c><see cref="Node::past_"/> == nullptr</c>) et <c>Node::trackedIndex_ &gt; 0</c>, déroule la chaîne via <c>future_</c>
	/// et écrit : (Track, Plane, Id, X, Y, Z, Intensité, Surface). Met à jour <see cref="nTraj_"/> au passage.\n
	/// Effets de bord : écrit dans le tampon <see cref="tracks_"/>, met à jour <see cref="nTraj_"/>.\n
	/// Complexité : \f$\mathcal{O}(M)\f$ avec M le nombre de nœuds exportés.
	/// </remarks>
	size_t printConnectionsASCII_trc();

	/// <summary>Assigne un indice de couleur aux structures suivies, constant le long d'une même trajectoire.</summary>
	/// <remarks>
	/// Les nœuds tête (<c><see cref="Node::past_"/> == nullptr</c>) reçoivent un nouvel indice incrémental, propagé aux descendants via <see cref="Node::past_"/>.\n
	/// Effets de bord : modifie <c>Node::colorIndex_</c> pour les nœuds suivis, remet à zéro pour les non suivis.\n
	/// Complexité : \f$\mathcal{O}(N)\f$ sur l'ensemble des nœuds.
	/// </remarks>
	void setColorToTrackedStructures();

	/// <summary>Étiquette les trajectoires valides (durée ≥ <see cref="minLife_"/>) avec un identifiant de suivi <c>Node::trackedIndex_</c>.</summary>
	/// <remarks>
	/// Pour chaque nœud tête, compte la longueur de la chaîne via <c>Node::future_</c> ; si \f$\text{longueur} \ge \text{minLife}\f$
	/// alors l'ensemble de la chaîne reçoit un identifiant incrémental, sinon <c>Node::trackedIndex_</c> est laissé à 0.\n
	/// Effets de bord : modifie <c>Node::trackedIndex_</c> (et donc l'ensemble des nœuds d'une trajectoire).\n
	/// Complexité : \f$\mathcal{O}(T)\f$ avec T le nombre de chaînes.
	/// </remarks>
	void setIndexToTrackedStructures();
};
} // namespace palm_tracer::tracking
