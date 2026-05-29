///-------------------------------------------------------------------------------------------------
///
/// \file Tracks.hpp
/// \brief Définition de la classe <see cref="Tracks"/> permettant de manipuler et de reconnecter des trajectoires issues du suivi de points (tracking) en microscopie.
/// La classe fournit des méthodes de chargement, de recomposition et d'export des trajectoires au format tableau 2D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 30/09/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <cstdint>
#include <vector>

#include "ComputesParams.hpp"
#include "Point.hpp"
#include "TrackingPoint.hpp"

namespace palm_tracer::tracking {
/// <summary>Classe manipulant un ensemble de trajectoires.
/// Chaque trajectoire est représentée comme une suite de points (<see cref="Point"/>).
/// La classe permet de charger des trajectoires, de les reconnecter en cas de scintillement, et de les exporter dans un format tabulaire générique.
/// </summary>
class Tracks {
	using track_t = std::vector<Point>; ///< Type représentant une trajectoire individuelle (suite de points).

public:
	Tracks()  = default; ///< Constructeur par défaut de la classe <see cref="Tracks"/>.
	~Tracks() = default; ///< Destructeur par défaut de la classe <see cref="Tracks"/>.

	/// <summary>Charge un ensemble de trajectoires précalculé depuis un tableau 2D brut.</summary>
	/// <param name="input">Pointeur vers un tableau contigu de doubles.
	/// Chaque ligne correspond à un point et doit contenir <see cref="N_TRACKING"/> colonnes (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <param name="nRow">Nombre de lignes (points) dans le tableau.</param>
	/// <param name="pixelSize">Taille d'un pixel en unités réelles (nanomètre).</param>
	/// <remarks>Cette fonction initialise le vecteur interne de trajectoires <c>t_</c>.</remarks>
	void LoadTracking(const double* input, size_t nRow, double pixelSize);

	/// <summary>Exporte les trajectoires actuelles vers un tableau 2D brut.</summary>
	/// <param name="output">Pointeur vers un tableau de sortie contigu de doubles.
	/// Chaque ligne produite correspond à un point avec les mêmes colonnes que l'entrée (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <returns>Nombre total de valeurs écrites (soit 8 × nombre de points exportés).</returns>
	/// <remarks>La colonne “Track” correspond à l'index de la trajectoire dans le tableau interne, réindexé à partir de 1 après éventuelles opérations de reconnexion.</remarks>
	size_t FillOutput(double* output) const;

	/// <summary></summary>
	/// <param name="input"></param>
	/// <param name="output"></param>
	/// <param name="nPoints"></param>
	/// <param name="maxDistance"></param>
	/// <returns></returns>
	static uint64_t Tracking(const double* input, uint64_t* output, size_t nPoints, double maxDistance);

	/// <summary>Reconnecte des trajectoires disjointes lorsqu'elles sont séparées par un scintillement (frames manquantes).</summary>
	/// <param name="mode">Mode de dispersion utilisé pour calculer la distance maximale de recherche (<see cref="blinking_mode"/>).</param>
	/// <param name="maxDuration">Durée maximale (en nombre de plans) entre la fin d'une trajectoire et le début d'une autre pour qu'une reconnexion soit possible.</param>
	/// <param name="maxSpeed">Vitesse maximale (en pixels/plan) supposée pour un point entre deux plans consécutifs.</param>
	/// <remarks>
	/// La reconnexion est effectuée de manière gloutonne :
	/// - priorité au plus petit décalage temporel \f$\Delta t\f$,
	/// - puis, à  \f$\Delta t\f$ égal, priorité à la distance spatiale minimale,
	/// sous réserve que cette distance soit inférieure au seuil calculé par <c>getMaxDistance</c>.
	/// </remarks>
	void BlinkingReconnection(size_t mode, size_t maxDuration, double maxSpeed);

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
	[[nodiscard]] static std::vector<double> ComputeMSD(const track_t& track, size_t start = 0, size_t end = 0);

	/// <summary>Calcule le <b>MSD</b> pour l'ensemble des trajectoires.</summary>
	/// <returns>Tableau 2D avec un nombre de lignes égal au nombre de trajectoires ; chaque ligne a une taille égale au lag maximal de la trajectoire (les lags sans observations valent -1).</returns>
	/// <remarks>En cas de blinking, les trajectoires plus longues contiennent des trous. Le résultat est moins précis (bruité), mais c'est convenu.</remarks>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllMSD() const;

	/// <summary>Estime la <b>diffusion instantanée</b> (Instant D) le long d'une trajectoire, via des <b>fenêtres glissantes</b> de longueur <c>ComputesParams::fitLength</c>.
	/// Pour chaque fenêtre, on calcule le vecteur MSD local pour les décalages \f$i = 1..\text{fitLength}\f$, puis on ajuste une droite : \f$\mathrm{MSD}(\tau) \approx 2\,d\,D_{\text{inst}}\,\tau + C\f$, où \f$d\f$ est la dimension (2 ou 3) et \f$\tau = i\,\Delta t\f$.
	/// </summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Un vecteur de taille \f$\max(0, N - \text{fitLength})\f$ où \f$N\f$ est la longueur de la trajectoire.</returns>
	[[nodiscard]] static std::vector<double> ComputeInstantDiffusion(const track_t& track, const ComputesParams& p);

	/// <summary>Estime la <b>diffusion instantanée</b> pour l'ensemble des trajectoires.</summary>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Tableau 2D (une ligne par trajectoire) ; chaque ligne a une taille \f$\max(0,\,N-\text{fitLength})\f$ spécifique à la trajectoire.</returns>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllInstantDiffusion(const ComputesParams& p) const;

	/// <summary>Calcule des métriques de diffusion à partir d'une trajectoire et réalise un ou plusieurs fittings du MSD.</summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Tableau contenant le résultat des métriques (taille et contenu dépendant du mode d'ajustement).</returns>
	[[nodiscard]] static std::vector<double> ComputeFitMetrics(const track_t& track, const ComputesParams& p);

	/// <summary>Calcule des métriques de diffusion pour les trajectoires et réalise un ou plusieurs fittings du MSD.</summary>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Tableau 2D des métriques (une ligne par trajectoire), avec un format dépendant du <see cref="fit_mode"/> choisi.</returns>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllFitMetrics(const ComputesParams& p) const;

	/// <summary>Surcharge de l'opérateur de flux de sortie pour afficher les trajectoires.</summary>
	/// <param name="os">Flux de sortie.</param>
	/// <param name="tracks">Objet <see cref="Tracks"/> à afficher.</param>
	/// <returns>Référence vers le flux de sortie, pour chaînage.</returns>
	friend std::ostream& operator<<(std::ostream& os, const Tracks& tracks);

protected:
	std::vector<track_t> t_;                              ///< Ensemble des trajectoires manipulées.
	static constexpr double INSTANT_DIFFUSION_MIN = 1e-6; ///< Diffusion instantanée minimale.

	/// <summary></summary>
	/// <param name="points"></param>
	/// <param name="nPoints"></param>
	/// <returns></returns>
	[[nodiscard]] static std::vector<TrackingPoint> loadPoints(const double* points, size_t nPoints);

	/// <summary></summary>
	/// <param name="points"></param>
	/// <returns></returns>
	[[nodiscard]] static std::vector<PlaneRange> getPlaneRanges(const std::vector<TrackingPoint>& points);

	/// <summary></summary>
	/// <param name="points"></param>
	/// <param name="planes"></param>
	/// <param name="maxDistance"></param>
	static void track(std::vector<TrackingPoint>& points, const std::vector<PlaneRange>& planes, double maxDistance);

	/// <summary></summary>
	/// <param name="points"></param>
	/// <param name="output"></param>
	/// <returns></returns>
	static uint64_t fillTrackingIdOutput(const std::vector<TrackingPoint>& points, uint64_t* output);

	/// <summary> Calcule la distance maximale autorisée pour une reconnexion en fonction du <paramref name="mode"/> choisi, de la durée du scintillement et de la vitesse maximale spécifiée.</summary>
	/// <param name="mode">Mode de dispersion des points (<see cref="blinking_mode"/>).</param>
	/// <param name="duration">Durée du scintillement (\f$\Delta t\f$, en nombre de plans).</param>
	/// <param name="maxSpeed">Vitesse maximale supposée (micromètre/plan).</param>
	/// <returns>La distance maximale de recherche (en micromètre) pour ce couple (<paramref name="mode"/>, \f$\Delta t\f$).</returns>
	static double getMaxDistance(size_t mode, double duration, double maxSpeed);

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
	static size_t getFitResSize(fit_mode mode);
};
} // namespace palm_tracer::tracking
