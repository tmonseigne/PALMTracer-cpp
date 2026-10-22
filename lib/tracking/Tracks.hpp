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
#include <iosfwd>
#include <vector>

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

	/// <summary>Donne accès aux trajectoires stockées pour les traitements en place de la bibliothèque de tracking.</summary>
	/// <returns>Référence modifiable au conteneur interne, sans copie des trajectoires.</returns>
	[[nodiscard]] std::vector<Track>& GetTracks() { return t_; }

	/// <summary>Donne accès en lecture seule aux trajectoires stockées pour les analyses.</summary>
	/// <returns>Référence constante au conteneur interne, sans copie des trajectoires.</returns>
	[[nodiscard]] const std::vector<Track>& GetTracks() const { return t_; }

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

	/// <summary>Surcharge de l'opérateur de flux de sortie pour afficher les trajectoires.</summary>
	/// <param name="os">Flux de sortie.</param>
	/// <param name="tracks">Objet <see cref="Tracks"/> à afficher.</param>
	/// <returns>Référence vers le flux de sortie, pour chaînage.</returns>
	friend std::ostream& operator<<(std::ostream& os, const Tracks& tracks);

protected:
	std::vector<Track> t_; ///< Ensemble des trajectoires manipulées.

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
};
} // namespace palm_tracer::tracking
