///-------------------------------------------------------------------------------------------------
///
/// \file Tracks.hpp
/// <summary>Définition de la classe <see cref="Tracks"/> permettant de manipuler et de reconnecter des trajectoires issues du suivi de points (tracking) en
/// microscopie.
/// La classe fournit des méthodes de chargement, de recomposition et d'export des trajectoires au format tableau 2D.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 30/09/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <cstdint>
#include <utility>
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
public:
	Tracks()  = default; ///< Constructeur par défaut de la classe <see cref="Tracks"/>.
	~Tracks() = default; ///< Destructeur par défaut de la classe <see cref="Tracks"/>.

	/// <summary>Charge un ensemble de trajectoires précalculé depuis un tableau 2D brut.</summary>
	/// <param name="input">Pointeur vers un tableau contigu de doubles.
	/// Chaque ligne correspond à un point et doit contenir <see cref="tracking_columns::Count"/> colonnes (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <param name="nRow">Nombre de lignes (points) dans le tableau.</param>
	/// <param name="pixelSize">Taille d'un pixel en unités réelles (nanomètre).</param>
	/// <remarks>Cette fonction initialise le vecteur interne de trajectoires <c>t_</c>.</remarks>
	void LoadTracks(const double* input, size_t nRow, double pixelSize);

	/// <summary>Exporte les trajectoires actuelles vers un tableau 2D brut.</summary>
	/// <param name="output">Pointeur vers un tableau de sortie contigu de doubles.
	/// Chaque ligne produite correspond à un point avec les mêmes colonnes que l'entrée (Track, Plan, Id, X, Y, Z, Intensité, Surface).</param>
	/// <returns>Nombre total de valeurs écrites (soit 8 × nombre de points exportés).</returns>
	/// <remarks>La colonne “Track” correspond à l'index de la trajectoire dans le tableau interne, réindexé à partir de 1 après éventuelles opérations de
	/// reconnexion.</remarks>
	size_t FillOutput(double* output) const;

	/// <summary>Construit des trajectoires 2D à partir d'un tableau contigu de localisations simples, puis exporte le couple <c>(Id, Track)</c> de chaque
	/// point.</summary>
	/// <param name="input">Tableau d'entrée contigu de doubles contenant 4 valeurs par point : <c>Id</c>, <c>Plane</c>, <c>X</c>, <c>Y</c>.</param>
	/// <param name="output">Tableau de sortie contigu préalloué en entiers 64 bits. Il reçoit deux valeurs par point : <c>Id</c>, puis identifiant de
	/// trajectoire.</param>
	/// <param name="nPoints">Nombre de points présents dans <paramref name="input"/>.</param>
	/// <param name="maxDistance">Distance maximale autorisée entre deux points de deux plans consécutifs pour créer une liaison.</param>
	/// <returns>Nombre total de valeurs écrites dans <paramref name="output"/>, soit <c>2 × nPoints</c>.</returns>
	/// <remarks>
	/// Les numéros de trajectoires sont recalculés à partir de 1 après l'association gloutonne.
	/// L'entrée peut être non triée : les points sont triés en interne par plan puis par identifiant avant la recherche des liaisons.
	/// La fonction est pensée pour l'interface DLL/ctypes et ne conserve aucun état dans l'instance <see cref="Tracks"/>.
	/// </remarks>
	static uint64_t Tracking(const double* input, uint64_t* output, size_t nPoints, double maxDistance);

	/// <summary>Reconnecte des trajectoires disjointes lorsqu'elles sont séparées par un scintillement (frames manquantes).</summary>
	/// <param name="mode">Mode de dispersion utilisé pour calculer la distance maximale de recherche (<see cref="blinking_mode"/>).</param>
	/// <param name="maxDuration">Durée maximale (en nombre de plans) entre la fin d'une trajectoire et le début d'une autre pour qu'une reconnexion soit
	/// possible.</param>
	/// <param name="maxSpeed">Vitesse maximale (en pixels/plan) supposée pour un point entre deux plans consécutifs.</param>
	/// <remarks>
	/// La reconnexion est effectuée de manière gloutonne :
	/// - priorité au plus petit décalage temporel \f$\Delta t\f$,
	/// - puis, à  \f$\Delta t\f$ égal, priorité à la distance spatiale minimale,
	/// sous réserve que cette distance soit inférieure au seuil calculé par <c>getMaxDistance</c>.
	/// Un index temporel limite la génération parallèle des candidats à la fenêtre admissible.
	/// Chaque liste est triée de façon stable par getReconnectionCost, puis les connexions sont décidées séquentiellement
	/// dans l'ordre des têtes, en vérifiant les segments déjà pris. Les coûts NaN conservent le parcours historique.
	/// Dans un même plan, les indices d'origine sont parcourus dans l'ordre croissant pour préserver les ex aequo.
	/// Le classement des candidats admissibles est centralise dans <see cref="getReconnectionCost"/>.
	/// A coût exactement égal, le premier candidat rencontre est conservé ; l'association globale reste gloutonne.
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
	[[nodiscard]] static std::vector<double> ComputeMSD(const Track& track, size_t start = 0, size_t end = 0);

	/// <summary>Calcule le <b>MSD</b> pour l'ensemble des trajectoires.</summary>
	/// <returns>Tableau 2D avec un nombre de lignes égal au nombre de trajectoires ; chaque ligne a une taille égale au lag maximal de la trajectoire (les lags
	/// sans observations valent -1).</returns>
	/// <remarks>En cas de blinking, les trajectoires plus longues contiennent des trous. Le résultat est moins précis (bruité), mais c'est convenu.</remarks>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllMSD() const;

	/// <summary>Estime la <b>diffusion instantanée</b> (Instant D) le long d'une trajectoire, via des <b>fenêtres glissantes</b> de longueur
	/// <c>ComputesParams::fitLength</c>.
	/// Pour chaque fenêtre, on calcule le vecteur MSD local pour les décalages \f$i = 1..\text{fitLength}\f$, puis on ajuste une droite : \f$\mathrm{MSD}(\tau)
	/// \approx 2\,d\,D_{\text{inst}}\,\tau + C\f$, où \f$d\f$ est la dimension (2 ou 3) et \f$\tau = i\,\Delta t\f$.
	/// </summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Un vecteur de taille \f$\max(0, N - \text{fitLength})\f$ où \f$N\f$ est la longueur de la trajectoire.</returns>
	[[nodiscard]] static std::vector<double> ComputeInstantDiffusion(const Track& track, const ComputesParams& p);

	/// <summary>Estime la <b>diffusion instantanée</b> pour l'ensemble des trajectoires.</summary>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Tableau 2D (une ligne par trajectoire) ; chaque ligne a une taille \f$\max(0,\,N-\text{fitLength})\f$ spécifique à la trajectoire.</returns>
	[[nodiscard]] std::vector<std::vector<double>> ComputeAllInstantDiffusion(const ComputesParams& p) const;

	/// <summary>Calcule des métriques de diffusion à partir d'une trajectoire et réalise un ou plusieurs fittings du MSD.</summary>
	/// <param name="track">Trajectoire à analyser.</param>
	/// <param name="p">Paramètres de configuration des calculs.</param>
	/// <returns>Tableau contenant le résultat des métriques (taille et contenu dépendant du mode d'ajustement).</returns>
	[[nodiscard]] static std::vector<double> ComputeFitMetrics(const Track& track, const ComputesParams& p);

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
	std::vector<Track> t_;                                ///< Ensemble des trajectoires manipulées.
	static constexpr double INSTANT_DIFFUSION_MIN = 1e-6; ///< Diffusion instantanée minimale.

	/// <summary>Coût lexicographique d'une reconnexion admissible : frames minimales, puis distance minimale.</summary>
	/// <param name="frameGap">Différence de numéro de frame entre les deux extrémités, strictement positive.</param>
	/// <param name="distance">Distance spatiale entre les deux extrémités, dans l'unité de <c>Point::pos_</c>.</param>
	/// <returns>Paire \f$(\Delta f,d)\f$ a minimiser lexicographiquement ; ce n'est pas une somme pondérée.</returns>
	/// <remarks>Les limites de durée et de distance restent des contraintes d'admissibilité séparées du coût.</remarks>
	[[nodiscard]] static std::pair<size_t, double> getReconnectionCost(size_t frameGap, double distance);

	/// <summary>Candidat de reconnexion avec indice d'origine et coût déjà calculé.</summary>
	struct ReconnectionCandidate {
		size_t target;                  ///< Indice du segment cible dans t_.
		std::pair<size_t, double> cost; ///< Coût fourni uniquement par getReconnectionCost.
	};

	/// <summary>Génère et trie les reconnexions admissibles d'un segment, indépendamment des conflits gloutons.</summary>
	/// <param name="starts">Indices des segments non vides, triés par plan de début puis indice d'origine.</param>
	/// <param name="current">Indice du segment dont la fin est à prolonger.</param>
	/// <param name="mode">Mode de calcul de la distance maximale admissible.</param>
	/// <param name="maxDuration">Délai maximal inclusif, en plans.</param>
	/// <param name="maxSpeed">Paramètre spatial transmis à getMaxDistance.</param>
	/// <param name="candidates">Liste vidée puis remplie ; tous les candidats admissibles sont conservés.</param>
	/// <returns>Vrai si la liste est triée ; faux si un coût NaN impose de conserver l'ordre temporel initial.</returns>
	/// <remarks>Le tri stable utilise exclusivement getReconnectionCost. Les filtres taken et indice après la tête
	/// sont appliqués lors des connexions séquentielles. Les appels concurrents utilisent des listes de sortie distinctes.</remarks>
	[[nodiscard]] bool getReconnectionCandidates(const std::vector<size_t>& starts, size_t current, size_t mode,
	                                             size_t maxDuration, double maxSpeed, std::vector<ReconnectionCandidate>& candidates) const;

	/// <summary>Convertit un tableau brut <c>(Id, Plane, X, Y)</c> en liste de points de tracking triée par plan puis par identifiant.</summary>
	/// <param name="points">Tableau d'entrée contigu de doubles contenant 4 valeurs par point : <c>Id</c>, <c>Plane</c>, <c>X</c>, <c>Y</c>.</param>
	/// <param name="nPoints">Nombre de points à lire dans <paramref name="points"/>.</param>
	/// <returns>Vecteur de <see cref="TrackingPoint"/> initialisés, triés par <c>point.plane_</c> croissant puis par <c>point.id_</c> croissant.</returns>
	/// <remarks>
	/// Cette méthode ne conserve que les informations utiles au tracking 2D rapide : identifiant, plan et position XY en pixels.
	/// Les champs de chaînage (<c>previous</c>, <c>next</c>, <c>trackId</c>) conservent leurs valeurs par défaut définies par <see cref="TrackingPoint"/>.
	/// Complexité : \f$\mathcal{O}(N \log N)\f$ à cause du tri final.
	/// </remarks>
	[[nodiscard]] static std::vector<TrackingPoint> loadPoints(const double* points, size_t nPoints);

	/// <summary>Calcule les intervalles d'indices correspondant à chaque plan dans une liste de points déjà triée.</summary>
	/// <param name="points">Points triés par plan croissant puis par identifiant croissant.</param>
	/// <returns>
	/// Vecteur de plages <see cref="PlaneRange"/> indexé relativement au premier plan observé.
	/// Pour un plan <c>p</c>, la plage est accessible par <c>res[p - MinPlane]</c> et représente l'intervalle demi-ouvert <c>[begin, end[</c> dans <paramref
	/// name="points"/>.
	/// </returns>
	/// <remarks>
	/// Les plans absents sont représentés par une plage vide (<c>begin == end</c>) positionnée à l'indice d'insertion logique.
	/// Précondition : <paramref name="points"/> doit être trié par <c>point.plane_</c> croissant.
	/// Complexité : \f$\mathcal{O}(N + P)\f$, avec \f$N\f$ le nombre de points et \f$P\f$ le nombre de plans couverts entre le premier et le dernier plan.
	/// </remarks>
	[[nodiscard]] static std::vector<PlaneRange> getPlaneRanges(const std::vector<TrackingPoint>& points);

	/// <summary>Génère tous les candidats admissibles entre deux plans, sans tri ni décision d'association.</summary>
	/// <param name="points">Points triés par plan et identifiant.</param>
	/// <param name="current">Plage des points sources.</param>
	/// <param name="next">Plage des points cibles du plan suivant.</param>
	/// <param name="maxDistance">Rayon dont le carré définit le filtre exact historique.</param>
	/// <param name="candidates">Buffer vidé puis rempli avec tous les candidats admissibles, sans tri ; sa capacité est conservée.</param>
	/// <remarks>La recherche examine exhaustivement toutes les paires des deux plans.</remarks>
	static void getCandidates(const std::vector<TrackingPoint>& points, PlaneRange current, PlaneRange next,
	                          double maxDistance, std::vector<Candidate>& candidates);

	/// <summary>Associe gloutonnement les points de plans consécutifs pour construire des trajectoires 2D.</summary>
	/// <param name="points">Points triés et modifiables. Les champs <c>previous</c>, <c>next</c> et <c>trackId</c> sont mis à jour en place.</param>
	/// <param name="planes">Plages d'indices par plan, typiquement calculées par <see cref="getPlaneRanges"/>.</param>
	/// <param name="maxDistance">Distance maximale autorisée entre deux points de plans consécutifs.</param>
	/// <remarks>
	/// Pour chaque paire de plans consécutifs, tous les candidats dont la distance quadratique est inférieure à <c>maxDistance²</c> sont générés,
	/// triés par distance croissante puis par indices, puis validés de manière gloutonne si les deux extrémités sont encore libres.
	/// Les couples de plans sont indépendants pour la génération des candidats et peuvent être traités en parallèle, mais l'attribution finale des identifiants
	/// de trajectoire reste séquentielle pour garder un résultat déterministe.
	/// La recherche des candidats est exhaustive. L'indexation spatiale 2D évaluée sur la pile uTub complète
	/// (5,19 millions de localisations, jusqu'à 600 par frame) n'a pas démontré de gain significatif au regard
	/// des seuils pratiques et de la dispersion des mesures ; le parcours exhaustif est conservé pour sa simplicité.
	/// Chaque thread réutilise son vecteur de candidats entre transitions. La génération est quadratique par paire de plans.
	/// </remarks>
	static void track(std::vector<TrackingPoint>& points, const std::vector<PlaneRange>& planes, double maxDistance);

	/// <summary>Écrit les identifiants de trajectoire calculés dans un tableau de sortie compatible DLL.</summary>
	/// <param name="points">Points de tracking après association, contenant les champs <c>point.id_</c> et <c>trackId</c> finalisés.</param>
	/// <param name="output">Tableau de sortie contigu préalloué. Il reçoit deux entiers par point : <c>Id</c>, puis <c>Track</c>.</param>
	/// <returns>Nombre total de valeurs écrites dans <paramref name="output"/>, soit <c>2 × points.size()</c>.</returns>
	/// <remarks>
	/// L'ordre d'écriture suit l'ordre interne de <paramref name="points"/>, donc le tri par plan/identifiant effectué en amont.
	/// Le réappariement côté Python doit donc se faire par la colonne <c>Id</c> plutôt que par l'ordre des lignes d'origine.
	/// </remarks>
	static uint64_t fillTrackingIdOutput(const std::vector<TrackingPoint>& points, uint64_t* output);

	/// <summary> Calcule la distance maximale autorisée pour une reconnexion en fonction du <paramref name="mode"/> choisi, de la durée du scintillement et de
	/// la vitesse maximale spécifiée.</summary>
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
