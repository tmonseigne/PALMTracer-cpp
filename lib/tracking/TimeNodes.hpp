///-------------------------------------------------------------------------------------------------
///
/// \file TimeNodes.hpp
/// \brief Définition de la classe <see cref="TimeNodes"/>, conteneur des <see cref="Node"/> observés à un instant donné et
/// outils associés pour établir les voisinages temporels (\f$ t \leftrightarrow t + 1\f$) et les statistiques globales.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 18/04/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once

#include <vector>
#include "Node.hpp"
#include "Parameters.hpp"

namespace palm_tracer::tracking {
/**
 * <summary>
 *  Conteneur de nœuds  <see cref="Node"/> pour un instant temporel \f$t\f$, fournissant :
 *  - la création/gestion des <see cref="Node"/> du temps \f$t\f$ ;
 *  - l'établissement des voisinages avant/après (liaisons\f$ t \rightarrow t + 1\f$ et \f$ t \leftarrow t - 1\f$) sous une contrainte de vitesse ;
 *  - le tri des voisinages par intensité (priorisation heuristique) ;
 *  - l'agrégation de statistiques globales (sigma constants) ;
 *  - des opérations utilitaires d'initialisation et de perturbation aléatoire contrôlée.
 *
 *  Les unités sont héritées de <see cref="Node"/> et de <see cref="Parameters"/> :
 *  - <c><see cref="Node::Distance"/>(...)</c> retourne une distance spatiale en unités réelles cohérentes (p.ex. µm) ;
 *  - <c><see cref="Parameters::maxSpeed_"/></c> est une vitesse maximale (par pas de temps) dans la même unité que la distance.
 * </summary>
 * <remarks>
 *  <b>Complexités :</b>
 *  - Construction : \f$\mathcal{O}(nNodes)\f$ ;
 *  - <c><see cref="LinkNeighborsBeforeAndAfter"/></c> : \f$\mathcal{O}(nNodes \times nNodes_\mathrm{fut})\f$ ;
 *  - <c><see cref="SortByIntensityAllNeighbors"/></c> : \f$\mathcal{O}\!\left(\sum_i (d_i \times \log \times d_i)\right)\f$, avec \f$d_i\f$ la taille du voisinage du nœud \f$i\f$.
 * </remarks>
 */
class TimeNodes {
protected:
	size_t time_ = 0; ///< Indice temporel associé aux nœuds contenus (frame/plan) : \f$t\f$.

public:
	Parameters* params_ = nullptr; ///< Pointeur vers les paramètres globaux (non possédés).
	std::vector<Node> nodes_;      ///< Nœuds observés au temps \f$t\f$.
	size_t nNodes_ = 0;            ///< Nombre attendu/allocable de nœuds en \f$t\f$.

	/// <summary>Construit <c>nNodes</c> nœuds au temps <paramref name="time"/> et les réserve dans <see cref="nodes_"/>.</summary>
	/// <param name="time">Indice temporel \f$t\f$ (frame/plan) associé aux nœuds créés.</param>
	/// <param name="nNodes">Nombre de nœuds à instancier pour \f$t\f$.</param>
	/// <param name="params">Paramètres globaux (non possédés).</param>
	/// <remarks>Chaque <see cref="Node"/> est construit via <c><see cref="Node"/>(<paramref name="time"/>, <paramref name="params"/>)</c>.</remarks>
	TimeNodes(size_t time, size_t nNodes, Parameters* params);

	/// <summary>Destructeur (vide). Le conteneur <see cref="nodes_"/> est détruit automatiquement.</summary>
	~TimeNodes() { nodes_.clear(); }

	/// <summary>Accès direct au <paramref name="i"/>-ème nœud (sans vérification de borne).</summary>
	/// <param name="i">Indice dans <c>[0, nNodes_[</c>.</param>
	/// <returns>Référence modifiable vers le <see cref="Node"/> demandé.</returns>
	Node& operator[](const size_t i) { return nodes_[i]; }

	/// <summary>Établit les voisinages bidirectionnels entre chaque nœud présent (\f$t\f$) et les nœuds futurs (\f$t + 1\f$)
	/// dont la distance est strictement inférieure à <c>params_->maxSpeed_</c>.</summary>
	/// <param name="futures">Référence sur la structure <see cref="TimeNodes"/> du temps \f$t + 1\f$.</param>
	/// <remarks>
	///  Pour un couple \f$(p, f)\f$ (présent, futur), si \f$ \mathrm{Distance}(p, f) < \texttt{maxSpeed} \f$ alors :
	///  <c>p.neighbors_[1]</c> reçoit <c>&f</c> et <c>f.neighbors_[0]</c> reçoit <c>&p</c>.
	///  Les tableaux <c>neighbors_[0]</c> et <c>neighbors_[1]</c> portent respectivement les liens \f$ t - 1 \to t\f$ et \f$t \to t + 1\f$.
	/// </remarks>
	void LinkNeighborsBeforeAndAfter(TimeNodes& futures);

	/// <summary>Trie, pour chaque nœud, l'ensemble de ses voisinages \f$(t - 1)\f$ et \f$(t + 1)\f$ par intensité décroissante (<see cref="Node::intensity_"/>).</summary>
	void SortByIntensityAllNeighbors();

	/// <summary>Agrège, sur l'ensemble des nœuds du temps \f$t\f$, des statistiques nécessaires au calcul de sigmas constants (intensité, moyenne, distance).</summary>
	/// <param name="intensity">Accumulateur de la somme d'intensités.</param>
	/// <param name="nInt">Compteur associé à <paramref name="intensity"/>.</param>
	/// <param name="mean">Accumulateur de moyenne (au sens de <see cref="Node::ComputeConstantSigmas"/>).</param>
	/// <param name="nMean">Compteur associé à <paramref name="mean"/>.</param>
	/// <param name="distance">Accumulateur de distances (au sens de <see cref="Node::ComputeConstantSigmas"/>).</param>
	/// <param name="nDist">Compteur associé à <paramref name="distance"/>.</param>
	/// <remarks>
	///  Cette méthode se contente d'agréger les contributions élémentaires de chaque <see cref="Node"/>
	///  via <c><see cref="Node::ComputeConstantSigmas"/>(...)</c>.
	/// </remarks>
	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const;

	/// <summary>Lance, pour chaque nœud, l'initialisation de la connexion (heuristique de départ).</summary>
	void MakeAllInitializationConnections() { for (auto& node : nodes_) { node.MakeInitializationConnection(); } }

	/// <summary>Propose un changement aléatoire local sur un nœud, choisi uniformément dans <c>[0, nNodes_[</c>, et retourne la variation de coût associée (au sens de <see cref="Node::MakeChange"/>).</summary>
	/// <param name="limit">Ampleur maximale de la perturbation.</param>
	/// <returns>Variation de coût retournée par <c><see cref="Node::MakeChange"/>(<paramref name="limit"/>)</c>. Vaut <c>0.0</c> si <c>nNodes_ == 0</c>.</returns>
	double TryRandomChange(const double limit) { return nNodes_ == 0 ? 0.0 : nodes_[params_->RandomModulo(nNodes_)].MakeChange(limit); }
};
} // namespace palm_tracer::tracking
