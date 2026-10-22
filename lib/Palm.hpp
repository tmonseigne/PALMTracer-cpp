///-------------------------------------------------------------------------------------------------
///
/// \file Palm.hpp
/// <summary>Classe principale pour la détection de points lumineux dans une image 2D en niveaux de gris,
/// basée sur une décomposition en ondelettes et un post-traitement par segmentation.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 26/03/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <cstdint>

#include "Constants.hpp"
#include "Structs.hpp"

namespace palm_tracer {
/// <summary>Classe principale pour la détection de points lumineux dans une image 2D en niveaux de gris,
/// basée sur une décomposition en ondelettes et un post-traitement par segmentation.</summary>
///
/// <remarks> La classe Palm effectue les étapes suivantes :
///		- Décomposition en ondelettes de type ATrous.
///		- Seuillage adaptatif du plan sélectionné.
///		- Watershed pour la séparation des objets connectés.
///		- Calcul de caractéristiques (intensité, volume, centroides).
///		- Filtrage morphologique basé sur les seuils d'intensité et de volume.
///		- Fitting gaussien (optionnel) pour affiner la position des objets.
///
/// Le résultat est un tableau de points, chacun contenant <see cref="localization_columns::Count"/> doubles représentant
/// des informations telles que la position, l'intensité, les paramètres de l'ajustement, etc.
///
/// Le traitement est lancé via la méthode <see cref="Processing()"/>.</remarks>
class Palm {
public:
	/// <summary>Constructeur de la classe <see cref="Palm"/>.</summary>
	///
	/// Initialise les structures internes pour effectuer une détection basée sur une décomposition en ondelettes et une segmentation intelligente.
	///
	/// <param name="image">Pile contiguë de pixels non signés 16 bits, de taille X × Y × Z, ordonnée par plan.</param>
	/// <param name="points">Pointeur vers le tableau de sortie où seront stockés les points détectés (x, y, intensité, etc.).</param>
	/// <param name="pointsListSize">Taille totale allouée pour le tableau de points (doit être un multiple de <see cref="localization_columns::Count"/>).</param>
	/// <param name="dim">Dimensions de la pile (X, Y, Z).</param>
	/// <param name="thresholdVal">Seuil de détection appliqué sur les coefficients du plan d'ondelettes choisi.</param>
	/// <param name="watershedRatio">Seuil de fusion des bassins : 0 conserve les maxima distincts, 10 autorise les fusions.</param>
	/// <param name="fit">Mode d'ajustement à appliquer :
	///		- 0 : aucun.
	///		- 1 : Gaussien 2D avec sigma fixe.
	///		- 2 : Gaussien 2D avec sigma variable.
	///		- 3 : Gaussien 3D avec theta fixe.
	///		- 4 : Gaussien 3D complet.</param>
	/// <param name="fitParams"> paramètres de l'ajustement
	/// En cas d'ajustement Gaussien :
	///		- Taille de la ROI (fenêtre carrée centrée) utilisé pour l'ajustement autour de chaque point détecté.
	///		- Valeur initiale de sigma X, en pixels.
	///		- Valeur initiale de sigma Y, en pixels.
	///		- Valeur initiale de l'orientation theta, en radians.
	///	En cas d'ajustement par Splines :
	///		- Taille de la ROI (fenêtre carrée centrée) utilisé pour l'ajustement autour de chaque point détecté.
	///		- Pas sur Z en nanomètre.
	///		- Position initiale du Z en nanomètre.
	///		- Tableau de Coefficients d'une PSF optimale ? pour chaque Z ? (tableau 4D TODO Voir quelles dimensions est pour quoi)
	///		- Image de fond pour le CMOS et 0 pour le EMCD (pour définir le bruit standard). (taille standard d'une image)
	/// </param>
	Palm(const uint16_t* image, double* points, size_t pointsListSize, const Vec3I& dim,
	     double thresholdVal, double watershedRatio, size_t fit, const double* fitParams);

	/// <summary>Destructeur de la classe <see cref="Palm"/>.</summary>
	~Palm() = default;

	/// <summary>Lance le traitement complet de détection sur l'image d'entrée.
	///
	/// Étapes principales :
	/// - Décomposition en ondelettes (ATrous).
	/// - Seuillage du plan choisi et création d'une carte binaire de segmentation.
	/// - Application d'un Watershed intelligent.
	/// - Calcul des caractéristiques des objets segmentés.
	/// - Application éventuelle de filtres morphologiques (volume, intensité).
	/// - Fitting gaussien optionnel autour des objets détectés.
	/// - Remplissage du tableau <see cref="points_"/> avec les résultats (coordonnées, intensité, sigma, etc.).
	///
	/// Le tableau est rempli par blocs de <see cref="localization_columns::Count"/> doubles, sans ligne sentinelle.
	/// La capacité est répartie également entre les plans ; les détections excédentaires d'un plan sont ignorées.
	/// Chaque plan exécute son pipeline complet : ondelettes, segmentation, fitting puis tri Y/X.
	/// </summary>
	/// <returns>Le nombre de doubles écrits, soit le nombre de détections retenues multiplié par <see cref="localization_columns::Count"/>.</returns>
	[[nodiscard]] int Processing() const;

	/// <summary>Ajuste des positions fournies sans détection ni segmentation, en conservant leur ordre.</summary>
	/// <param name="image">Pile contiguë de pixels 16 bits, ordonnée par plan.</param>
	/// <param name="positions">Tableau contigu de <c>nPoints × 6</c> doubles : identifiant global, plan entier à partir de 1, X et Y en pixels à partir de 0,
	/// intensité de segmentation et surface de segmentation en pixels carrés.</param>
	/// <param name="nPoints">Nombre de positions à ajuster.</param>
	/// <param name="output">Tampon distinct des entrées, au format habituel de 18 colonnes par localisation.</param>
	/// <param name="outputSize">Capacité du tampon en doubles ; doit couvrir toutes les positions.</param>
	/// <param name="dim">Dimensions de la pile (X, Y, Z), strictement positives.</param>
	/// <param name="fit">Mode d'ajustement historique : 1 à 4 pour les gaussiennes, 5 pour les splines EMCCD.</param>
	/// <param name="fitParams">Paramètres habituels du modèle choisi ; ROI positive et compatible avec l'image.</param>
	/// <returns>Nombre de doubles écrits ; 0 si l'entrée est vide.</returns>
	/// <remarks>Les entrées respectent les contrats du pipeline de localisation. Les identifiants globaux sont recopiés depuis l'entrée ;
	/// les indices locaux suivent l'ordre d'apparition dans chaque plan. L'intensité et la surface de segmentation sont recopiées depuis l'entrée.
	/// Une ROI débordant de l'image conserve les coordonnées initiales et reçoit les sentinelles habituelles d'échec du fit.
	/// La validité des positions, des dimensions, des tampons et des paramètres du modèle est à la charge de l'appelant.</remarks>
	[[nodiscard]] static size_t FitLocalizations(const uint16_t* image, const double* positions, size_t nPoints,
	                                           double* output, size_t outputSize, const Vec3I& dim, size_t fit, const double* fitParams);

	/// <summary>Calcule automatiquement un seuil adapté à l'image en fonction des localisations détectées.
	///
	/// La méthode procède par itérations successives afin d'affiner le seuil de segmentation.
	/// Elle utilise un algorithme basé sur la détection de points lumineux (localisations)
	/// et la caractérisation statistique des pixels hors des régions d'intérêt (ROI) définies autour de ces localisations.
	///
	/// Étapes principales :
	/// - Initialisation de l'écart type (seuil) à partir de l'image brute.
	/// - À chaque itération :
	///   - Détection des points lumineux sur l'image avec le seuil courant via le processeur de frame interne.
	///   - Construction d'un masque binaire qui couvre les zones autour des localisations détectées (ROI).
	///   - Extraction des pixels hors masque (hors ROI).
	///   - Calcul de l'écart type des pixels hors ROI.
	///   - Mise à jour du seuil si l'écart type a changé significativement.
	/// - Arrêt si l'écart type ne varie plus ou si les pixels hors ROI sont insuffisants.
	///
	/// Le masque binaire est construit avec des zones carrées centrées sur chaque localisation, de taille définie par `fitParams_[0]` (taille ROI).
	/// </summary>
	/// <returns>Le seuil calculé correspondant à l'écart type final des pixels hors régions segmentées.</returns>
	[[nodiscard]] double AutoThreshold();

private:
	// Basic members
	const uint16_t* stack_;     ///< Image 3D d'entrée (plan x hauteur x largeur) en niveaux de gris.
	double* points_;            ///< Tableau de points de sortie (rempli avec les résultats de détection).
	Vec3I dim_;                 ///< Dimensions de la pile (X, Y, Z).
	size_t size_,               ///< Taille totale d'une image (height × width).
	       nPointsMaxPerPlane_; ///< Nombre maximal de localisations par plan, selon la capacité totale et le nombre de plans.

	// Décomposition members
	double thresholdVal_,   ///< Seuil de détection utilisé après la décomposition en ondelettes.
	       watershedRatio_; ///< Rapport utilisé dans l'algorithme de Watershed.

	size_t fit_;              ///< Mode de fitting à appliquer (0 = aucun, 1 à 4 = variantes 2D ou 3D de l'ajustement Gaussien, 5 = Ajustement par Splines).
	const double* fitParams_; ///< Paramètres de l'ajustement
};
} // namespace palm_tracer
