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
#include <cstdint>
#include <vector>

#include "Point.hpp"
#include "TracksFit.hpp"

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
	using track_t = std::vector<Point>; ///< Type représentant une trajectoire individuelle (suite de points).

public:
	Tracks()  = default; ///< Constructeur par défaut de la classe <c>Tracks</c>.
	~Tracks() = default; ///< Destructeur par défaut de la classe <c>Tracks</c>.

	/// <summary>Charge un ensemble de trajectoires précalculé depuis un tableau 2D brut.</summary>
	/// <param name="input">Pointeur vers un tableau contigu de doubles.
	/// Chaque ligne correspond à un point et doit contenir <see cref="N_TRACKING"/> colonnes (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <param name="nRow">Nombre de lignes (points) dans le tableau.</param>
	/// <param name="pixelSize">Taille d'un pixel en unitées réelles (nanomètre).</param>
	/// <remarks>Cette fonction initialise le vecteur interne de trajectoires <c>t_</c>.</remarks>
	void LoadTracking(const double* input, size_t nRow, double pixelSize);

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

	/// <summary>>Calcule le **Mean Squared Displacement (MSD)**, en considérant tous les décalages temporels (écart de plan) possibles sur une trajectoire.
	/// - Soit une trajectoire de longueur N : positions successives P₁, P₂, …, Pₙ.
	/// - Pour chaque décalage (ou *lag*) i = 1 … (plan(Pₙ) - plan(P₁)),
	///	  on calcule la moyenne des distances quadratiques : MSD(i) = mean_j (‖ P_{j+i} - P_j ‖²)
	///   avec j parcourant toutes les paires séparées de i pas.
	/// - La dernière colonne (lag = N-1) est donc basée sur un seul terme (‖Pₙ - P₁‖²).
	///
	/// En cas de Blinking, il est possible qu'un lag ne possède aucun ensemble de point le représentant, le MSD sera à -1.
	/// </summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="start">Premier point de la trajectoire à prendre (utile lorsqu'un sous ensemble de points est sélectionné pour le fit).</param>
	/// <param name="end">Dernier point de la trajectoire à prendre (utile lorsqu'un sous ensemble de points est sélectionné pour le fit).</param>
	/// <returns>Un vecteur de longueur \f$L_i=\mathrm{plane}_{\max}-\mathrm{plane}_{\min}\f$ (ou vide si la trajectoire contient moins de 2 points)</returns>
	/// <remarks>
	/// - Hypothèse : La trajectoire est triée par plan strictement croissant (aucun doublon de plan).
	/// - Unités : \c Point::GetDistance() fixe l’unité (nm). Appliquer un facteur d’échelle en amont si nécessaire.
	/// - Les bacs (lags) sans paires observées sont mis à <c>-1.0</c>.
	/// - Complexité : \f$\mathcal{O}\!\left(M^2\right)\f$ pour \f$M\f$ le nombre de points par trajectoire.
	/// </remarks>
	[[nodiscard]] static std::vector<double> ComputeMSD(const track_t& track, size_t start = 0, size_t end = 0);

	/// <summary>Calcule le **Mean Squared Displacement (MSD)**, en considérant tous les décalages temporels possibles.
	/// - Soit une trajectoire de longueur N : positions successives P₁, P₂, …, Pₙ.
	/// - Pour chaque décalage (ou *lag*) i = 1 … N-1,
	///   on calcule la moyenne des distances quadratiques : MSD(i) = mean_j (‖ P_{j+i} - P_j ‖²)
	///   avec j parcourant toutes les paires séparées de i pas.
	/// - La dernière colonne (lag = N-1) est donc basée sur un seul terme (‖Pₙ - P₁‖²).
	///
	/// </summary>
	/// <returns>Tableau 2D avec un nombre de lignes égal au nombre de trajectoires et chaque ligne de taille (N-1)
	/// (avec N la différence de plan maximale de la trajectoire, potentiellement supérieur au nombre de points de la trajectoire en cas de blinking).</returns>
	/// <remarks>En cas de blinking, les trajectoires plus longues contiennent des trous. Le résultat est moins précis (bruité), mais c'est convenu.</remarks>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllMSD() const;

	/// <summary>Estime la <b>diffusion instantanée</b> (Instant D) le long de la trajectoire, via des <b>fenêtres glissantes</b> de longueur <c>FitLength</c>.
	///
	///	Pour chaque fenêtre, on calcule le vecteur MSD local pour les décalages
	/// i = 1 … FitLength (dans la fenêtre), puis on ajuste une droite : MSD(τ) ≈ 2·d·D_inst·τ + C
	/// où d = dimension (2 ou 3), τ = i·Δt. La pente donne D_inst = slope / (2·d).
	/// </summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="fitLength">Taille de la fenêtre glissante (sur le MSD) pour le calcul de la pente.</param>
	/// <param name="is3D"> Si <c>true</c>, utilise d = 3 (MSD en 3D avec x,y,z) ; sinon d = 2 (x,y).</param>
	/// <returns></returns>
	[[nodiscard]] static std::vector<double> ComputeInstantDiffusion(const track_t& track, size_t fitLength, bool is3D);


	/// <summary>Estime la <b>diffusion instantanée</b> (Instant D) le long des trajectoires, via des <b>fenêtres glissantes</b> de longueur <c>FitLength</c>.\n
	/// Pour chaque fenêtre, on calcule le vecteur MSD local pour les décalages
	/// i = 1 … FitLength (dans la fenêtre), puis on ajuste une droite : MSD(τ) ≈ 2·d·D_inst·τ + C
	/// où d = dimension (2 ou 3), τ = i·Δt. La pente donne D_inst = slope / (2·d).</summary>
	/// <param name="fitLength">Taille de la fenêtre glissante pour le calcul de la pente.</param>
	/// <param name="is3D"> Si <c>true</c>, utilise d = 3 (MSD en 3D avec x,y,z) ; sinon d = 2 (x,y).</param>
	/// <returns>Tableau 2D avec un nombre de lignes égal au nombre de trajectoires et chaque ligne de taille (max(0, 1 + floor((N_t - (fit_length+1)) / stride))) (N la longueur de la trajectoire).</returns>
	/// <remarks> Pour chaque fenêtre, on prend les FitLength premiers points du MSD local (lags 1..FitLength).
	/// On les régréssionne vs τ = (1..FitLength)·Δt, puis on renvoie D_inst = slope / (2·d), avec d = 2 (2D) ou 3 (3D).</remarks>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllInstantDiffusion(size_t fitLength, bool is3D) const;


	/// <summary>Calcule des métriques de diffusion à partir d'une trajectoire et réalise un ou plusieurs fittings du MSD.</summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="is3D"> Si <c>true</c>, utilise d = 3 (MSD en 3D avec x,y,z) ; sinon d = 2 (x,y).</param>
	/// <param name="mode">Le mode d'ajustement (voir <see cref="fit_mode"/>).</param>
	/// <param name="params">Paramètres de configuration de l'ajustement.</param>
	/// <returns>Tabelau contenant le résultat des métriques (taille et contenu dépendant du mode d'ajustement).</returns>
	[[nodiscard]] static std::vector<double> ComputeFitMetrics(const track_t& track, bool is3D, fit_mode mode, const std::vector<double>& params);


	/// <summary>Calcule des métriques de diffusion pour les trajectoires et réalise un ou plusieurs fittings du MSD.</summary>
	/// <param name="is3D"> Si <c>true</c>, utilise d = 3 (MSD en 3D avec x,y,z) ; sinon d = 2 (x,y).</param>
	/// <param name="mode">Le mode d'ajustement (voir <see cref="fit_mode"/>).</param>
	/// <param name="params">Paramètres de configuration de l'ajustement.</param>
	/// <returns>Tabelau contenant le résultat des métriques (taille et contenu dépendant du mode d'ajustement).</returns>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllFitMetrics(bool is3D, fit_mode mode, const std::vector<double>& params) const;

	/// <summary>Surcharge de l’opérateur de flux de sortie pour afficher les trajectoires.</summary>
	/// <param name="os">Flux de sortie.</param>
	/// <param name="tracks">Objet <c>Tracks</c> à afficher.</param>
	/// <returns>Référence vers le flux de sortie, pour chaînage.</returns>
	friend std::ostream& operator<<(std::ostream& os, const Tracks& tracks);

protected:
	std::vector<track_t> t_;                              ///< Ensemble des trajectoires manipulées.
	double pixelSize_                             = 160;  ///< Taille d'un pixel en unitées réelles (nanomètre).
	static constexpr double INSTANT_DIFFUSION_MIN = 1e-6; ///< Diffusion instantannée minimale.


	/// <summary> Calcule la distance maximale autorisée pour une reconnexion en fonction du mode choisi, de la durée du scintillement et de la vitesse maximale spécifiée.</summary>
	/// <param name="mode">Mode de dispersion des points (<see cref="blinking_mode"/>).</param>
	/// <param name="duration">Durée du scintillement (Δt, en nombre de plans).</param>
	/// <param name="maxSpeed">Vitesse maximale supposée (pixels/plan).</param>
	/// <returns> La distance maximale de recherche (en pixels) pour ce couple (mode, Δt). </returns>
	static double getMaxDistance(size_t mode, double duration, double maxSpeed);

	/// <summary>Retourne le nombre d'éléments dnas le tableau final attendu pour le <see cref="fit_mode"/> courant.</summary>
	/// <returns>
	/// \f[
	/// \begin{cases}
	/// 6 & \text{pour } \mathtt{fit\_mode::linear} \quad ([Da, Db, Derr, a, b, err])\\[2pt]
	/// 7 & \text{pour } \mathtt{fit\_mode::power} \quad ([Da, Db, Derr, α, a, err, \mu(\text{dist})])\\[2pt]
	/// 8 & \text{pour } \mathtt{fit\_mode::exponential} \quad ([Da, Db, Derr, a, \tau, b, err, \sqrt{\frac{3}{4}\times(b+a)}])
	/// \end{cases}
	/// \f]
	/// </returns>
	/// <remarks>En cas de  <c>mode_</c> non reconnu, retourne <c>0</c> (à éviter en pratique).</remarks>
	static size_t getFitResSize(fit_mode mode);
};
} // namespace palm_tracer
