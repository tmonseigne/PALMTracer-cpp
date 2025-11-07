///-------------------------------------------------------------------------------------------------
///
/// \file Parameters.hpp
/// \brief Paramètres globaux de configuration et utilitaires aléatoires pour le suivi/les mesures.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 18/04/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <random>

namespace palm_tracer {
/// <summary>
/// Conteneur de paramètres généraux utilisés par les modules de suivi et de calculs, ainsi que quelques utilitaires liés au générateur pseudo-aléatoire.
///
/// Points d'attention :
/// - Les grandeurs de vitesse sont exprimées par défaut en micromètre/plan (cf. conventions de <see cref="Tracking"/>).
/// - Les bornes spatiales <see cref="maxX_"/>, <see cref="maxY_"/>, <see cref="maxZ_"/> stockent des tailles (et non des indices), dans les mêmes unités que les positions manipulées par <see cref="Node"/>.
/// - Le générateur <see cref="rng_"/> est initialisé avec une graine fixe (42) pour assurer la reproductibilité des tests.</summary>
class Parameters {
public:
	Parameters()  = default; ///< Constructeur par défaut de la classe <see cref="Parameters"/>.
	~Parameters() = default; ///< Destructeur par défaut de la classe <see cref="Parameters"/>.

	// ==================== Constantes globales (configuration par défaut) ====================
	/// <summary>Distance maximale de recherche (par défaut) pour certaines heuristiques de liaison.
	/// Unité : micromètre (ou unité spatiale cohérente avec <see cref="Node"/>).</summary>
	static constexpr double DIST_MAX = 20.0;

	static constexpr size_t N_STEPS      = 5;   ///< Nombre d'étapes discrètes par itération pour des procédures incrémentales.
	static constexpr size_t N_ITERATIONS = 300; ///< Nombre maximal d'itérations autorisé pour une procédure d'optimisation/raffinement.
	static constexpr size_t MODEL        = 2;   ///< Identifiant du modèle sélectionné (convention interne).

	// ==================== Paramètres de coût / contraintes (suivi) ====================
	double maxSpeed_{};       ///< Vitesse maximale admissible lors de la prédiction/l'association (micromètre/plan).
	double costBirthDeath_{}; ///< Coût unitaire de naissance/mort d'une piste dans l'algorithme d'association (sans dimension).

	// ==================== Statistiques pré-calculées (normalisation des coûts) ====================
	double stdIntensity_{};     ///< Écart-type de l'intensité des détections (même unité que <see cref="Node::intensity_"/>).
	double stdMean_{};          ///< Écart-type de la moyenne d'intensité sur une fenêtre (unité de <see cref="Node::intensity_"/>).
	double stdDistance_{};      ///< Écart-type des distances entre points (même unité spatiale que <see cref="Node"/>).
	double covMeanIntensity_{}; ///< Covariance entre la moyenne d'intensité et l'intensité instantanée (unités combinées).

	// ==================== Bornes spatiales du domaine observé ====================
	double maxX_{}; ///< Taille/borne maximale suivant X (unité spatiale de <see cref="Node"/>).
	double maxY_{}; ///< Taille/borne maximale suivant Y (unité spatiale de <see cref="Node"/>).
	double maxZ_{}; ///< Taille/borne maximale suivant Z (même unité; p. ex. micromètre).

	// ==================== Générateur pseudo-aléatoire ====================
	std::mt19937 rng_{42}; ///< Générateur Mersenne-Twister initialisé avec une graine fixe (42) pour la reproductibilité.

	// ==================== Utilitaires ====================
	/// <summary>Retourne un entier uniforme dans \f$\{0,1,\dots,n-1\}\f$.</summary>
	/// <param name="n">Module strictement positif (taille de l'intervalle).</param>
	/// <returns>Une valeur entière dans l'intervalle \f$[0,n[\f$.</returns>
	/// <remarks>
	/// Précondition : \f$n>0\f$. La distribution est obtenue par réduction modulo de <see cref="rng_()"/>.
	/// Pour éviter le biais de modulo lorsqu'il est critique, préférer <c>std::uniform_int_distribution</c>.
	/// </remarks>
	[[nodiscard]] size_t RandomModulo(const size_t n) { return rng_() % n; }

	/// <summary>Retourne un flottant uniforme sur \f$[0,1]\f$.</summary>
	/// <returns>Une valeur dans \f$[0,1]\f$ issue de la normalisation du tirage de <see cref="rng_"/>.</returns>
	/// <remarks>Distribution continue \f$U(0,1)\f$ obtenue par division par <c>std::mt19937::max()</c>.</remarks>
	[[nodiscard]] double Randomfloat() { return static_cast<double>(rng_()) / static_cast<double>(std::mt19937::max()); }
};
} // namespace palm_tracer
