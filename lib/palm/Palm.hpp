///-------------------------------------------------------------------------------------------------
///
/// \file Palm.hpp
/// \brief  Classe principale pour la détection de points lumineux dans une image 2D en niveaux de gris,
/// basée sur une décomposition en ondelettes et un post-traitement par segmentation.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <cstdint>
#include <vector>

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
/// Le résultat est un tableau de points, chacun contenant N_LOCALIZATIONS doubles représentant
/// des informations telles que la position, l'intensité, les paramètres de l'ajustement, etc.
///
/// Le traitement est lancé via la méthode Processing().</remarks>
class Palm {
public:
	/// <summary>Constructeur de la classe <c>Palm</c>.</summary>
	///
	/// Initialise les structures internes pour effectuer une détection basée sur une décomposition en ondelettes et une segmentation intelligente.
	///
	/// <param name="image">Image d'entrée au format tableau 2D (grayscale, taille height × width).</param>
	/// <param name="points">Pointeur vers le tableau de sortie où seront stockés les points détectés (x, y, intensité, etc.).</param>
	/// <param name="pointsListSize">Taille totale allouée pour le tableau de points (doit être un multiple de N_LOCALIZATIONS).</param>
	/// <param name="dim">Dimensions de la pile (X, Y, Z).</param>
	/// <param name="thresholdVal">Seuil de détection appliqué sur les coefficients du plan d'ondelettes choisi.</param>
	/// <param name="watershedRatio">Ratio utilisé par l'algorithme de Watershed pour la séparation des objets.</param>
	/// <param name="fit">Mode d'ajustement à appliquer :
	///		- 0 : aucun.
	///		- 1 : Gaussien 2D avec sigma fixe.
	///		- 2 : Gaussien 2D avec sigma variable.
	///		- 3 : Gaussien 3D avec theta fixe.
	///		- 4 : Gaussien 3D complet.</param>
	/// <param name="fitParams"> paramètres de l'ajustement
	/// En cas d'ajustement Gaussien :
	///		- Taille du ROI (fenêtre carrée centrée) utilisé pour l'ajustement autour de chaque point détecté.
	///		- Valeur initiale de sigma utilisé dans l'ajustement (sigmaX = sigma et sigmaY = 2* sigma).
	///		- Valeur initiale de l'orientation theta pour l'ajustement.
	///	En cas d'ajustement par Splines :
	///		- Taille du ROI (fenêtre carrée centrée) utilisé pour l'ajustement autour de chaque point détecté.
	///		- Pas sur Z en nanomètre.
	///		- Position initiale du Z en nanomètre.
	///		- Tableau de Coefficients d'une PSF optimale ? pour chaque Z ? (tableau 4D TODO Voir quelles dimensions est pour quoi)
	///		- Image de fond pour le CMOS et 0 pour le EMCD (pour définir le bruit standard). (taille standard d'une image)
	/// </param>
	Palm(const uint16_t* image, double* points, size_t pointsListSize, const Vec3I& dim,
	     double thresholdVal, double watershedRatio, size_t fit, const double* fitParams);

	/// <summary>Destructeur de la classe <c>Palm</c>.</summary>
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
	/// - Remplissage du tableau `points_` avec les résultats (coordonnées, intensité, sigma, etc.).
	///
	/// Le tableau `points_` est rempli avec des blocs de N_LOCALIZATIONS double, représentant les paramètres pour chaque détection. Une ligne vide (valeurs = -1) est ajoutée à la fin.
	/// </summary>
	/// <returns>Le nombre total de valeurs (doubles) écrites dans `points_`, y compris la ligne finale vide. (i.e. `nb_points * N_LOCALIZATIONS + N_LOCALIZATIONS`) </returns>
	[[nodiscard]] int Processing() const;

	/// <summary>Calcule automatiquement un seuil adapté à l'image en fonction des localisations détectées.
	///
	/// La méthode procède par itérations successives afin d’affiner le seuil de segmentation.
	/// Elle utilise un algorithme basé sur la détection de points lumineux (localisations)
	/// et la caractérisation statistique des pixels hors des régions d’intérêt (ROI) définies autour de ces localisations.
	///
	/// Étapes principales :
	/// - Initialisation de l’écart type (seuil) à partir de l’image brute.
	/// - À chaque itération :
	///   - Détection des points lumineux sur l’image avec le seuil courant via `processingPlane`.
	///   - Construction d’un masque binaire qui couvre les zones autour des localisations détectées (ROI).
	///   - Extraction des pixels hors masque (hors ROI).
	///   - Calcul de l’écart type des pixels hors ROI.
	///   - Mise à jour du seuil si l’écart type a changé significativement.
	/// - Arrêt si l’écart type ne varie plus ou si les pixels hors ROI sont insuffisants.
	///
	/// Le masque binaire est construit avec des zones carrées centrées sur chaque localisation, de taille définie par `fitParams_[0]` (taille ROI).
	/// </summary>
	/// <returns>Le seuil calculé correspondant à l’écart type final des pixels hors régions segmentées.</returns>
	[[nodiscard]] double AutoThreshold();

private:
	// Basic members
	const uint16_t* stack_;     ///< Image 3D d'entrée (plan x hauteur x largeur) en niveaux de gris.
	double* points_;            ///< Tableau de points de sortie (rempli avec les résultats de détection).
	Vec3I dim_;                 ///< Dimensions de la pile (X, Y, Z).
	size_t size_,               ///< Taille totale d'une image (height × width).
	       nPointsMaxPerPlane_; ///< Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).

	// Decomposition members
	double thresholdVal_,   ///< Seuil de détection utilisé après la décomposition en ondelettes.
	       watershedRatio_; ///< Rapport utilisé dans l'algorithme de Watershed.

	size_t fit_;              ///< Mode de fitting à appliquer (0 = aucun, 1 à 4 = variantes 2D ou 3D de l'ajustement Gaussien, 5 = Ajustement par Splines).
	const double* fitParams_; ///< Paramètres de l'ajustement

	/// <summary>Lance le traitement complet de détection sur le plan passé en paramèter de l'image d'entrée.</summary>
	///
	/// Étapes principales :
	/// - Décomposition en ondelettes (ATrous).
	/// - Seuillage du plan choisi et création d'une carte binaire de segmentation.
	/// - Application d'un Watershed intelligent.
	/// - Calcul des caractéristiques des objets segmentés.
	/// - Application éventuelle de filtres morphologiques (volume, intensité).
	/// - Fitting gaussien optionnel autour des objets détectés.
	///
	/// Le tableau `points_` est rempli avec des blocs de N_LOCALIZATIONS double, représentant les paramètres pour chaque détection. Une ligne vide (valeurs = -1) est ajoutée à la fin.
	///
	/// <param name="plane">Index du plan à analyser.</param>
	/// <param name="points">Résultat de l'analyse.</param>
	void processingPlane(size_t plane, std::vector<std::array<double, N_LOCALIZATIONS>>& points) const;
};
} // namespace palm_tracer
