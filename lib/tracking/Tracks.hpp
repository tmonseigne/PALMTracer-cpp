///-------------------------------------------------------------------------------------------------
///
/// \file Tracks.hpp
/// \brief Définition de la classe \c Tracks permettant de manipuler et de reconnecter des trajectoires issues du suivi de points (tracking) en microscopie.
/// La classe fournit des méthodes de chargement, de recomposition et d’export des trajectoires au format tableau 2D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <vector>

#include "Point.hpp"

namespace palm_tracer {
/// <summary>Modes de croissance de la distance maximale de recherche lors de la reconnexion de trajectoires séparées par scintillement.</summary>
enum class blinking_mode : uint8_t {
	immobile,  ///< Distance maximale constante avec le temps (recherche dans un rayon fixe).
	diffusive, ///< Distance maximale croissante en racine carrée du temps écoulé.
	linear     ///< Distance maximale croissante linéairement avec le temps écoulé.
};

/// <summary>Classe manipulant un ensemble de trajectoires.
/// Chaque trajectoire est représentée comme une suite de points (<see cref="Point"/>).
/// La classe permet de charger des trajectoires, de les reconnecter en cas de scintillement, et de les exporter dans un format tabulaire générique.
/// </summary>
class Tracks {
public:
	Tracks()  = default; ///< Constructeur par défaut de la classe <c>Tracks</c>.
	~Tracks() = default; ///< Destructeur par défaut de la classe <c>Tracks</c>.

	/// <summary>Charge un ensemble de trajectoires précalculé depuis un tableau 2D brut.</summary>
	/// <param name="input">Pointeur vers un tableau contigu de doubles.
	/// Chaque ligne correspond à un point et doit contenir <see cref="N_TRACKING"/> colonnes (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <param name="nRow">Nombre de lignes (points) dans le tableau.</param>
	/// <remarks>Cette fonction initialise le vecteur interne de trajectoires <c>t_</c>.</remarks>
	void LoadTracking(const double* input, size_t nRow);

	/// <summary>Exporte les trajectoires actuelles vers un tableau 2D brut.</summary>
	/// <param name="output">Pointeur vers un tableau de sortie contigu de doubles.
	/// Chaque ligne produite correspond à un point avec les mêmes colonnes que l’entrée (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <returns>Nombre total de valeurs écrites (soit 8 × nombre de points exportés).</returns>
	/// <remarks>La colonne “Track” correspond à l’index de la trajectoire dans le tableau interne, réindexé à partir de 1 après éventuelles opérations de reconnexion.</remarks>
	size_t FillOutput(double* output) const;

	// TODO À terme, implémenter une méthode de suivi plus avancée que l'approche aléatoire.
	// void Tracking(const double* points, double* tracks, const size_t minLife, const double decrease, const size_t ntimes, Parameters* params);

	/// <summary>Reconnecte des trajectoires disjointes lorsqu’elles sont séparées par un scintillement (frames manquantes).</summary>
	/// <param name="mode">Mode de dispersion utilisé pour calculer la distance maximale de recherche (<see cref="blinking_mode"/>).</param>
	/// <param name="maxDuration">Durée maximale (en nombre de plans) entre la fin d’une trajectoire et le début d’une autre pour qu’une reconnexion soit possible.</param>
	/// <param name="maxSpeed">Vitesse maximale (en pixels/plan) supposée pour un point entre deux plans consécutifs.</param>
	/// <remarks>
	/// La reconnexion est effectuée de manière gloutonne :
	/// - priorité au plus petit décalage temporel Δt,
	/// - puis, à Δt égal, priorité à la distance spatiale minimale,
	/// sous réserve que cette distance soit inférieure au seuil calculé par <c>getMaxDistance</c>.
	/// </remarks>
	void BlinkingReconnection(size_t mode, size_t maxDuration, double maxSpeed);

	/// <summary>Surcharge de l’opérateur de flux de sortie pour afficher les trajectoires.</summary>
	/// <param name="os">Flux de sortie.</param>
	/// <param name="tracks">Objet <c>Tracks</c> à afficher.</param>
	/// <returns>Référence vers le flux de sortie, pour chaînage.</returns>
	friend std::ostream& operator<<(std::ostream& os, const Tracks& tracks);

protected:
	using track_t = std::vector<Point>; ///< Type représentant une trajectoire individuelle (suite de points).
	std::vector<track_t> t_;            ///< Ensemble des trajectoires manipulées.


	/// <summary> Calcule la distance maximale autorisée pour une reconnexion en fonction du mode choisi, de la durée du scintillement et de la vitesse maximale spécifiée.</summary>
	/// <param name="mode">Mode de dispersion des points (<see cref="blinking_mode"/>).</param>
	/// <param name="duration">Durée du scintillement (Δt, en nombre de plans).</param>
	/// <param name="maxSpeed">Vitesse maximale supposée (pixels/plan).</param>
	/// <returns> La distance maximale de recherche (en pixels) pour ce couple (mode, Δt). </returns>
	static double getMaxDistance(size_t mode, double duration, double maxSpeed);
};
} // namespace palm_tracer
