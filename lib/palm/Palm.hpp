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
#include <vector>

#include "Constants.hpp"

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
/// Le résultat est un tableau de points, chacun contenant N_SEGMENTS doubles représentant
/// des informations telles que la position, l'intensité, les paramètres de l'ajustement, etc.
///
/// Le traitement est lancé via la méthode Processing(). </remarks>
class Palm {
public:
	/// <summary>Constructeur de la classe <c>Palm</c>.</summary>
	///
	/// Initialise les structures internes pour effectuer une détection basée sur une décomposition en ondelettes et une segmentation intelligente.
	///
	/// <param name="image">Image d'entrée au format tableau 2D (grayscale, taille height × width).</param>
	/// <param name="points">Pointeur vers le tableau de sortie où seront stockés les points détectés (x, y, intensité, etc.).</param>
	/// <param name="pointsListSize">Taille totale allouée pour le tableau de points (doit être un multiple de N_SEGMENTS).</param>
	/// <param name="height">Hauteur d'une image (en pixels).</param>
	/// <param name="width">Largeur d'une image (en pixels).</param>
	/// <param name="planes">Nombre de plans de la pile</param>
	/// <param name="thresholdVal">Seuil de détection appliqué sur les coefficients du plan d'ondelettes choisi.</param>
	/// <param name="watershedRatio">Ratio utilisé par l'algorithme de Watershed pour la séparation des objets.</param>
	/// <param name="fit">Mode d'ajustement à appliquer :
	///		- 0 : aucun.
	///		- 1 : Gaussien 2D avec sigma fixe.
	///		- 2 : Gaussien 2D avec sigma variable.
	///		- 3 : Gaussien 3D avec theta fixe.
	///		- 4 : Gaussien 3D complet.</param>
	/// <param name="sigma">Valeur initiale de sigma utilisée dans le fitting (sigmaX = sigma et sigmaY = 2* sigma).</param>
	/// <param name="theta">Valeur initiale de l'orientation theta pour le fitting.</param>
	/// <param name="roi">Taille du ROI (fenêtre carrée centrée) utilisé pour le fitting autour de chaque point détecté.</param>
	Palm(const unsigned short* image, double* points, const size_t pointsListSize, const size_t height, const size_t width, const size_t planes,
	     const double thresholdVal, const double watershedRatio, const size_t fit, const double sigma, const double theta, const size_t roi)
		: image_(image), points_(points), height_(height), width_(width), planes_(planes), size_(height * width),
		  nPointsMax_(pointsListSize / N_SEGMENTS), nPointsMaxPerPlane_((pointsListSize / N_SEGMENTS) / planes),
		  thresholdVal_(thresholdVal), watershedRatio_(watershedRatio), fit_(fit), roi_(roi), initSigmaX_(sigma), initSigmaY_(2 * sigma), initTheta_(theta) { }

	/// <summary>Destructeur de la classe <c>Palm</c>.</summary>
	~Palm() = default;

	/// <summary>Lance le traitement complet de détection sur l'image d'entrée.</summary>
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
	/// Le tableau `points_` est rempli avec des blocs de N_SEGMENTS double, représentant les paramètres pour chaque détection. Une ligne vide (valeurs = -1) est ajoutée à la fin.
	///
	/// <returns>Le nombre total de valeurs (doubles) écrites dans `points_`, y compris la ligne finale vide. (i.e. `nb_points * N_SEGMENTS + N_SEGMENTS`) </returns>
	int Processing() const;

private:
	// Basic members
	const unsigned short* image_; ///< Image 3D d'entrée (plan x hauteur x largeur) en niveaux de gris.
	double* points_;              ///< Tableau de points de sortie (rempli avec les résultats de détection).
	size_t height_,               ///< Hauteur de l'image.
	       width_,                ///< Largeur de l'image.
	       planes_,               ///< Nombre de plans de la pile.
	       size_,                 ///< Taille totale d'une image (height × width).
	       nPointsMax_,           ///< Nombre maximum de points (pour éviter le débordement dans le tableau `points_`).
	       nPointsMaxPerPlane_;   ///< Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).

	// Decomposition members
	double thresholdVal_,   ///< Seuil de détection utilisé après la décomposition en ondelettes.
	       watershedRatio_; ///< Rapport utilisé dans l'algorithme de Watershed.

	// Gaussian fit members
	size_t fit_,        ///< Mode de fitting à appliquer (0 = aucun, 1 à 4 = variantes 2D ou 3D de l'ajustement Gaussien).
	       roi_;        ///< Taille du ROI utilisé pour le fitting gaussien (fenêtre autour du point).
	double initSigmaX_, ///< Valeur initiale de sigmaX utilisée dans le fitting gaussien.
	       initSigmaY_, ///< Valeur initiale de sigmaY utilisée dans le fitting gaussien.
	       initTheta_;  ///< Valeur initiale de theta (orientation) utilisée dans le fitting gaussien.

	// Splines members


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
	/// Le tableau `points_` est rempli avec des blocs de N_SEGMENTS double, représentant les paramètres pour chaque détection. Une ligne vide (valeurs = -1) est ajoutée à la fin.
	///
	/// <param name="plane">Index du plan à analyser.</param>
	/// <param name="points">Résultat de l'analyse.</param>
	void processingPlane(size_t plane, std::vector<std::array<double, N_SEGMENTS>>& points) const;
};
