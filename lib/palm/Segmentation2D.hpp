///-------------------------------------------------------------------------------------------------
///
/// \file Segmentation2D.hpp
/// \brief  Classe de traitement d'images segmentées en 2D, permettant l'allocation de structures,
/// le calcul de statistiques et le filtrage basé sur des critères morphologiques et d'intensité.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include "ATrous.hpp"
#include "SegmentationInfo.hpp"

namespace palm_tracer {
/// <summary>Classe de traitement d'images segmentées en 2D, permettant l'allocation de structures,
/// le calcul de statistiques et le filtrage basé sur des critères morphologiques et d'intensité. </summary>
///
/// Cette classe manipule des structures 2D identifiées dans une image segmentée et stockées sous forme de labels.
/// Elle permet notamment :
///		- L'allocation et la désallocation d'un vecteur de structures (`SegmentationInfo`).
///		- Le calcul de statistiques géométriques et photométriques par structure.
///		- Le filtrage de structures selon des critères définis par l'utilisateur.
///		- Le traitement de fusion locale par carte de potentiel (`filterStructureWS`).
///
/// Elle repose sur une image segmentée (labelisée), ainsi qu'éventuellement sur des cartes de données ou de potentiels.
class Segmentation2D {
public:
	/// <summary></summary>
	/// <param name="decomposition"></param>
	/// <param name="threshold">Seuil de détection utilisé après la décomposition en ondelettes.</param>
	/// <param name="ratio">Seuil relatif utilisé pour filtrer les germes faibles (si 0, le filtrage est ignoré).</param>
	Segmentation2D(const ATrous& decomposition, double threshold, double ratio);

	/// <summary></summary>
	~Segmentation2D() = default;

	/// <summary>Calcule les statistiques par structure à partir de l'image d'intensité et de la segmentation.</summary>
	void ComputeStatistics();

	/// <summary>Filtre les structures segmentées selon divers critères fournis en paramètres.
	///
	/// Le filtrage s'effectue en deux passes :
	/// 1. Critères morphologiques et photométriques (taille, intensité, rapport de forme…).
	/// 2. Suppression de structures proches de structures plus brillantes (si activé via `useParam[6]`).
	///
	/// L'image segmentée est mise à jour en conséquence (réindexation des structures conservées).</summary>
	void FilterData();

	/// <summary>Retourne la liste des structures segmentées avec leurs informations (centroïde, intensité, etc.).</summary>
	/// <returns>Un vecteur d'objets SegmentationInfo contenant les caractéristiques des structures détectées.</returns>
	/// <remarks>Le mot-clé [[nodiscard]] indique que le résultat ne doit pas être ignoré.</remarks>
	[[nodiscard]] std::vector<SegmentationInfo> Infos() const { return infos_; }

	/// <summary> Effectue une segmentation intelligente de type watershed.
	/// Cette méthode combine une détection des maxima locaux, un premier watershed, un filtrage morphologique optionnel, puis un second watershed.</summary>
	void WatershedIntelligent();

private:
	static constexpr double THRESHOLD = 1e-8; ///< Seuil minimal utilisé pour certaines comparaisons flottantes.
	size_t height_,                           ///< Hauteur de l'image d'entrée.
	       width_,                            ///< Largeur de l'image d'entrée.
	       size_;                             ///< Taille totale de l'image (height_ * width_).
	double ratio_;                            ///< Rapport utilisé dans l'algorithme de Watershed.
	std::vector<SegmentationInfo> infos_;     ///< Vecteur de structures 2D analysées ou filtrées.
	std::vector<uint16_t> map_;               ///< Image segmentée (carte binaire des objets détectés).
	std::vector<double> output_;              ///< Carte des valeurs utilisées pour la segmentation (après seuillage).

	/// <summary>Alloue un vecteur de structures (`infos_`) pour un nombre donné de segments.</summary>
	/// <param name="nbStruct">Nombre de structures initiales à allouer.</param>
	void alloc(size_t nbStruct);

	/// <summary>Applique un filtrage des structures par fusion locale via la carte de potentiel.
	///
	/// Cette fonction fusionne les structures en contact si la différence de potentiel justifie une absorption locale.
	///
	/// </summary>
	void filterStructureWS();

	/// <summary>Détecte les maxima locaux dans la carte de potentiel et initialise une carte de germes pour le watershed.
	/// Chaque maximum local reçoit un label unique, étendu à ses voisins immédiats s'ils ne sont pas dominés.</summary>
	void getLocalMax();

	/// <summary> Applique l'algorithme de segmentation watershed à partir des germes préalablement détectés.
	/// Utilise une file de priorité pour propager les labels des régions à travers la carte de potentiel. </summary>
	void watershed();
};
} // namespace palm_tracer
