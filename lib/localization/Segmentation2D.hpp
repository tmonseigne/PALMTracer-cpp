///-------------------------------------------------------------------------------------------------
///
/// \file Segmentation2D.hpp
/// \brief  Classe de traitement d'images segmentées en 2D, permettant l'allocation de structures,
/// le calcul de statistiques et le filtrage basé sur des critères morphologiques et d'intensité.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 26/03/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <vector>

#include "ATrous.hpp"
#include "SegmentationInfo.hpp"

namespace palm_tracer::localization {
/// <summary>Classe de traitement d'images segmentées en 2D, permettant l'allocation de structures,
/// le calcul de statistiques et le filtrage basé sur des critères morphologiques et d'intensité.</summary>
///
/// Cette classe manipule des structures 2D identifiées dans une image segmentée et stockées sous forme de labels.
/// Elle permet notamment :
///		- L'allocation et la désallocation d'un vecteur de structures (<see cref="SegmentationInfo"/>).
///		- Le calcul de statistiques géométriques et photométriques par structure.
///		- Le filtrage de structures selon des critères définis par l'utilisateur.
///		- Le traitement de fusion locale par carte de potentiel (<see cref="filterStructureWS()"/>).
///
/// Elle repose sur une image segmentée (labelisée), ainsi qu'éventuellement sur des cartes de données ou de potentiels.
class Segmentation2D {
public:
	/// <summary>Segmente le plan d'ondelettes seuillé, puis mesure et filtre les régions obtenues.</summary>
	/// <param name="decomposition">Décomposition dont le plan <c>PLANE_ID</c> fournit le potentiel de segmentation.</param>
	/// <param name="threshold">Seuil de détection utilisé après la décomposition en ondelettes.</param>
	/// <param name="ratio">Seuil relatif de fusion des bassins, et non un pourcentage.
	/// 0 conserve les germes distincts (Watershed activé dans Python) ; 10 autorise leur absorption par un voisin plus intense
	/// (Watershed désactivé dans Python). Une valeur élevée favorise donc les fusions, pas la séparation.</param>
	Segmentation2D(const ATrous& decomposition, double threshold, double ratio);

	/// <summary>Libère les cartes et les mesures des régions.</summary>
	~Segmentation2D() = default;

	/// <summary>Calcule les statistiques par structure à partir de l'image d'intensité et de la segmentation.</summary>
	void ComputeStatistics();

	/// <summary>Applique les filtres morphologiques et photométriques activés dans la configuration compilée.
	/// Les régions conservées et leurs étiquettes sont ensuite compactées dans la carte segmentée.</summary>
	/// <remarks>La configuration actuelle active uniquement la surface minimale : une région doit contenir
	/// strictement plus de <c>VOLUME_MIN</c> pixels. Les autres filtres restent désactivés.</remarks>
	void FilterData();

	/// <summary>Accède aux mesures des régions segmentées sans copier leur vecteur.</summary>
	/// <returns>Référence constante vers les informations détenues par la segmentation.</returns>
	/// <remarks>La segmentation doit rester vivante et ne pas être modifiée pendant l'utilisation de la référence.
	/// Une opération de segmentation, de filtrage ou de recalcul peut invalider les références sur ses éléments.
	/// Utiliser <c>auto snapshot = segmentation.Infos()</c> pour conserver une copie indépendante.</remarks>
	[[nodiscard]] const std::vector<SegmentationInfo>& Infos() const noexcept { return infos_; }

	/// <summary>Effectue une segmentation intelligente de type watershed.
	/// Cette méthode combine une détection des maxima locaux, un premier watershed, un filtrage morphologique optionnel, puis un second watershed.</summary>
	void WatershedIntelligent();

private:
	static constexpr double THRESHOLD = 1e-8; ///< Seuil minimal utilisé pour certaines comparaisons flottantes.
	Vec2I dim_;                               ///< Dimensions de l'image (X, Y).
	size_t size_;                             ///< Taille totale de l'image (height_ * width_).
	double ratio_;                            ///< Rapport utilisé dans l'algorithme de Watershed.
	std::vector<SegmentationInfo> infos_;     ///< Vecteur de structures 2D analysées ou filtrées.
	std::vector<uint32_t> map_;               ///< Labels des régions, sans limite à 65535 objets par image.
	std::vector<double> output_;              ///< Carte des valeurs utilisées pour la segmentation (après seuillage).

	/// <summary>Alloue un vecteur de structures (`infos_`) pour un nombre donné de segments.</summary>
	/// <param name="nbStruct">Nombre de structures initiales à allouer.</param>
	void alloc(size_t nbStruct);

	/// <summary>Fusionne les bassins en contact lorsque le contraste relatif de leur frontière est inférieur au seuil.</summary>
	/// <remarks>Pour un bassin de potentiels minimal m, maximal M et maximal en frontière B, le critère est
	/// \f$R=(M-B)/(M-m) &lt; \mathrm{ratio}\f$, avec un voisin strictement plus intense et une dynamique non nulle.
	/// Comme \f$0\leq R\leq 1\f$, le seuil historique 10 ne protège pas les vallées profondes.
	/// Les étiquettes sont positives ; 0 désigne l'absence de région, jamais l'indice de la première région.</remarks>
	void filterStructureWS();

	/// <summary>Détecte les maxima locaux dans la carte de potentiel et initialise une carte de germes pour le watershed.
	/// Chaque maximum local reçoit un label unique, étendu à ses voisins immédiats s'ils ne sont pas dominés.</summary>
	void getLocalMax();

	/// <summary>Applique l'algorithme de segmentation watershed à partir des germes préalablement détectés.
	/// Utilise une file de priorité pour propager les labels des régions à travers la carte de potentiel.</summary>
	void watershed();
};
} // namespace palm_tracer::localization
