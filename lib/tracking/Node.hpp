///-------------------------------------------------------------------------------------------------
///
/// \file Node.hpp
/// \brief Définition de la classe <see cref="Node"/> et API publique pour l'algorithme de réassignation temporelle.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 18/04/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <vector>

#include "Parameters.hpp"

namespace palm_tracer::tracking {
/// <summary>Représente une détection ponctuelle à l'instant <see cref="time_"/> et porte les liens temporels (<see cref="past_"/>, <see cref="future_"/>) vers les noeuds voisins.
/// La classe fournit :
/// - des mesures locales (position, intensité, moyenne, surface, distance aux bords…),
/// - des métriques de similarité (*likeness*) combinant distance, intensité, moyenne d'intensité,
/// - des opérations élémentaires de réassignation des liens (déconnexion, permutation, double permutation, etc.) sous critère de type "Metropolis".
/// </summary>
class Node {
public:
	//==================== Données observables ====================
	double x_{},         ///< Coordonnées X en unités réelles (plan image).
	       y_{},         ///< Coordonnées Y en unités réelles (plan image).
	       z_{},         ///< Coordonnées Z en unités réelles (plan image).
	       intensity_{}, ///< Intensité totale.
	       mean_{},      ///< Intensité moyenne.
	       surf_{},      ///< Surface projetée (aire, 2D).
	       id_{};        ///< Identifiant externe (optionnel).

	//==================== Contexte ====================
	Parameters* params_; ///< Contexte algorithme / paramètres globaux (non possédé).

	//==================== Chaînage temporel ====================
	size_t time_{};                      ///< Indice temporel discret (frame).
	Node* future_{};                     ///< Lien vers le nœud à \f$t + 1\f$ (ou <c>nullptr</c> si libre).
	Node* past_{};                       ///< Lien vers le nœud à \f$t - 1\f$ (ou <c>nullptr</c> si libre).
	size_t trackedIndex_{};              ///< Index d'appartenance à une trajectoire (si applicable).
	std::array<size_t, 2> colorIndex_{}; ///< Indices de coloration (UI / debug).

	// Voisinages par type (0: local, 1: candidats temporels)
	std::array<std::vector<Node*>, 2> neighbors_{};

	//==================== Pré-calculs (mis à jour via <see cref="UpdateConstants"/>) ====================
	double distanceToBorder_{},  ///< Distance minimale aux bords du FOV.
	       likeness_{},          ///< Similarité intrinsèque du nœud (utilisée quand <c>n == nullptr</c>).
	       likenessDistance_{},  ///< Contribution distance pré-calculée.
	       likenessIntensity_{}, ///< Contribution intensité pré-calculée.
	       likenessMean_{};      ///< Contribution intensité moyenne pré-calculée.

	//==================== Cycle de vie ====================
	/// <summary>Construit un nœud à l'instant <c>time</c> avec un pointeur vers <see cref="Parameters"/>.</summary>
	/// <param name="time">Indice temporel discret.</param>
	/// <param name="params">Paramètres globaux (non possédés).</param>
	Node(const size_t time, Parameters* params) : params_(params), time_(time) { }

	/// <summary>Destructeur : nettoie les voisinages et réinitialise l'état de traçage.</summary>
	~Node();

	//==================== Mesures / Similarités ====================
	/// <summary>Distance euclidienne 2D au nœud <c>node</c>.</summary>
	/// <param name="node">nœud cible.</param>
	/// <returns>\f$\sqrt{(x-x')^2+(y-y')^2}\f$.</returns>
	[[nodiscard]] double Distance(const Node& node) const;

	/// <summary>Calcule la similarité globale (*likeness*) avec un nœud <c>n</c>.
	/// Si <c>n == nullptr</c>, renvoie la similarité intrinsèque de ce nœud pénalisée par le coût naissance/mort.</summary>
	/// <param name="n">nœud cible ou <c>nullptr</c> pour la valeur locale.</param>
	/// <returns>Score de similarité adimensionnel.</returns>
	double Likeness(const Node* n = nullptr) const;

	/// <summary>Met à jour <see cref="likenessDistance_"/>, <see cref="likenessIntensity_"/>, <see cref="likenessMean_"/> et <see cref="likeness_"/> à partir de l'état courant.</summary>
	void UpdateConstants();

	/// <summary>Accumule des estimateurs de variabilité locale afin d'estimer des \f$\sigma\f$ constants
	/// pour intensité, moyenne et distance.</summary>
	/// <param name="intensity">Accumulateur de \f$\sum (\Delta I)^2\f$.</param>
	/// <param name="nInt">Compteur de paires pour l'intensité.</param>
	/// <param name="mean">Accumulateur de \f$\sum (\Delta m)^2\f$.</param>
	/// <param name="nMean">Compteur de paires pour la moyenne.</param>
	/// <param name="distance">Accumulateur de \f$\sum (\Delta d)^2\f$.</param>
	/// <param name="nDist">Compteur de paires pour la distance.</param>
	void ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const;

	/// <summary>Calcule la distance minimale aux bords du FOV et l'enregistre dans <see cref="distanceToBorder_"/>.</summary>
	void ComputeBorderDistance();

	//==================== Connexions (initialisation + MCMC) ====================
	/// <summary>Lie ce nœud à son meilleur candidat temporel (dans <c><see cref="neighbors_"/>[1]</c>) si disponible, en initialisation gloutonne.</summary>
	void MakeInitializationConnection();

	/// <summary>Propose et évalue une modification locale de graphe par "simple linkage".
	/// La variation d'énergie est acceptée suivant un critère de type <c>Metropolis</c>.</summary>
	/// <param name="limit">Température effective (seuil d'acceptation).</param>
	/// <returns>Variation d'énergie (\f$\Delta E\f$) si acceptée, ou \f$0\f$ sinon.</returns>
	double MakeChange(const double limit) { return makeChangeSimpleLinkage(limit); }

private:
	//==================== Outils internes ====================
	/// <summary>Test d'acceptation de type <b>Metropolis</b> pour une variation d'énergie \f$\Delta E\f$.</summary>
	/// <param name="deltaE">Variation d'énergie proposée \f$\Delta E\f$ (positive si la configuration se dégrade).</param>
	/// <param name="temperature">Température effective (contrôle la probabilité d'acceptation des dégradations).</param>
	/// <returns><c>true</c> si la proposition est acceptée, <c>false</c> sinon.</returns>
	[[nodiscard]] bool metropolis(double deltaE, double temperature) const;

	/// <summary>Rayon effectif 2D du nœud courant à partir de sa surface projetée.</summary>
	/// <returns>\f$r = \sqrt{\frac{S}{\pi}}\f$.</returns>
	/// <remarks> 2D : \f$r = \sqrt{\frac{S}{\pi}}\f$. Pour 3D, remplacer par \f$r = {\frac{3V}{4\pi}}^(1/3)\f$ si <c>GetVolume()</c> est disponible.</remarks>
	[[nodiscard]] double sphericEffect() const { return std::sqrt(surf_ / std::numbers::pi); }

	/// <summary>Rayon effectif combiné de deux nœuds (somme des rayons 2D).</summary>
	/// <param name="n1">Premier nœud.</param>
	/// <param name="n2">Second nœud.</param>
	/// <returns><c>n1.sphericEffect() + n2.sphericEffect()</c>.</returns>
	static double sphericEffect(const Node& n1, const Node& n2) { return n1.sphericEffect() + n2.sphericEffect(); }

	/// <summary>Distance de ce nœud au barycentre pondéré des intensités de <paramref name="n1"/> et <paramref name="n2"/>.</summary>
	/// <param name="n1">Premier nœud de référence.</param>
	/// <param name="n2">Second nœud de référence.</param>
	/// <returns>\f$\sqrt{(x - x_{12})^2 + (y - y_{12})^2}\f$, avec \f$(x_{12}, y_{12})\f$ le barycentre pondéré par l'intensité.</returns>
	[[nodiscard]] double distance(const Node& n1, const Node& n2) const;

	/// <summary>Retourne un voisin aléatoire parmi <see cref="neighbors_"/>[<paramref name="type"/>].</summary>
	/// <param name="type">Indice du groupe de voisinage (0 : local spatial, 1 : candidats temporels).</param>
	/// <returns>Pointeur vers un nœud voisin choisi uniformément.</returns>
	[[nodiscard]] Node* getRandomNeighbor(const size_t type = 0) const { return neighbors_[type][params_->RandomModulo(neighbors_[type].size())]; }

	//==================== Similarités élémentaires ====================
	/// <summary>Score élémentaire de similarité pour une différence scalaire.</summary>
	/// <param name="difference">Différence mesurée \f$\Delta\f$ (distance, intensité, moyenne...).</param>
	/// <param name="sigma">Échelle robuste \f$\sigma\f$ (écart-type ou proxy).</param>
	/// <returns>\f$1 - \dfrac{|\Delta|}{3\, \sigma}\f$.</returns>
	static double getLikeness(const double difference, const double sigma) { return 1. - std::abs(difference) / (3. * sigma); }


	/// <summary>Similarité adoucie en présence de clignotement (blinking).</summary>
	/// <param name="n">Nœud cible ou <c>nullptr</c> pour la valeur locale.</param>
	/// <returns>Score de similarité prenant en compte une interpolation avec les valeurs locales.</returns>
	double likenessBlink(const Node* n = nullptr) const;

	/// <summary>Terme distance de la similarité (contraste spatial).</summary>
	/// <param name="n">Nœud cible ou <c>nullptr</c> pour utiliser la distance au bord.</param>
	/// <returns>Score de similarité lié à la distance.</returns>
	double likenessDistance(const Node* n = nullptr) const;

	/// <summary>Terme intensité de la similarité.</summary>
	/// <param name="n">Nœud cible ou <c>nullptr</c> pour utiliser l'intensité locale.</param>
	/// <returns>Score de similarité lié à l'intensité.</returns>
	/// <remarks>Hypothèse : \f$\sigma_{I} \propto \frac{\sigma_{\text{noise}}}{\mu}\f$.</remarks>
	double likenessIntensity(const Node* n = nullptr) const;

	/// <summary>Terme moyenne d'intensité de la similarité.</summary>
	/// <param name="n">Nœud cible ou <c>nullptr</c> pour utiliser la moyenne locale.</param>
	/// <returns>Score de similarité lié à la moyenne d'intensité.</returns>
	/// <remarks>Hypothèse 2D : \f$\sigma_{\mu} \propto \text{surf}^{-1/6} \times \frac{\sigma_{\text{noise}}}{\mu}\f$.</remarks>
	double likenessMean(const Node* n = nullptr) const;

	//==================== Opérateurs de réassignation ("simple linkage") ====================
	/// <summary>Propose un changement local de connexions (déconnexion, permutation, etc.), évalué puis accepté via <see cref="metropolis"/>.</summary>
	/// <param name="limit">Température/limite d'acceptation pour le critère de <c>Metropolis</c>.</param>
	/// <returns>\f$\Delta E\f$ accepté (>0 amélioration, <0 dégradation acceptée), ou \f$0\f$ si refus.</returns>
	double makeChangeSimpleLinkage(double limit);

	/// <summary>Teste si un nœud candidat peut être relié sans conflit temporel.</summary>
	/// <param name="node">Candidat (même time que <see cref="future_"/>).</param>
	/// <returns><c>true</c> si non relié (<see cref="past_"/> nul) ou s'il s'agit déjà de <see cref="future_"/>.</returns>
	bool isSimpleLinkageAvailable(const Node* node) const { return node->past_ == nullptr || node == future_; }

	/// <summary>Proposition : couper le lien courant <see cref="future_"/> (déconnexion).</summary>
	/// <param name="limit">Température/limite d'acceptation.</param>
	/// <returns>Gain \f$\Delta E\f$ si accepté, \f$0\f$ sinon.</returns>
	[[nodiscard]] double simpleLinkageDisconnection(double limit) const;

	/// <summary>Proposition : réassigner <see cref="future_"/> vers un autre candidat libre.</summary>
	/// <param name="limit">Température/limite d'acceptation.</param>
	/// <param name="node">Nouveau candidat (libre).</param>
	/// <returns>Gain \f$\Delta E\f$ si accepté, \f$0\f$ sinon.</returns>
	double simpleLinkageReassignment(double limit, Node* node);

	/// <summary>Proposition : casser un conflit puis reconnecter en cascade.</summary>
	/// <param name="limit">Température/limite d'acceptation.</param>
	/// <param name="node">Candidat provoquant un conflit.</param>
	/// <returns>Gain \f$\Delta E\f$ si accepté, \f$0\f$ sinon.</returns>
	double simpleLinkageBreakAndReconnect(double limit, Node* node);

	/// <summary>Proposition : échanger les paires temporelles (swap à 2 nœuds).</summary>
	/// <param name="limit">Température/limite d'acceptation.</param>
	/// <param name="node">Candidat déjà relié.</param>
	/// <returns>Gain \f$\Delta E\f$ si accepté, \f$0\f$ sinon.</returns>
	double simpleLinkageSwap(double limit, Node* node);

	/// <summary>Proposition : double permutation avec reliaison d'un troisième nœud libre.</summary>
	/// <param name="limit">Température/limite d'acceptation.</param>
	/// <param name="node">Candidat initial.</param>
	/// <param name="count">Nombre de nœuds libres éligibles dans le voisinage du <see cref="past_"/> de <paramref name="node"/>.</param>
	/// <returns>Gain \f$\Delta E\f$ si accepté, \f$0\f$ sinon.</returns>
	double simpleLinkageDoubleSwap(double limit, Node* node, size_t count);
};
} // namespace palm_tracer::tracking
