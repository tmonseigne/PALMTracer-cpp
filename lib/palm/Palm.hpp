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

#include "ATrous.hpp"
#include "Constants.hpp"
#include "Segmentation2D.hpp"


/// <summary>Classe principale pour la détection de points lumineux dans une image 2D en niveaux de gris,
/// basée sur une décomposition en ondelettes et un post-traitement par segmentation.</summary>
///
/// <remarks>
/// La classe Palm effectue les étapes suivantes :
///		- Décomposition en ondelettes de type ATrous.
///		- Seuillage adaptatif du plan sélectionné.
///		- Watershed pour la séparation des objets connectés.
///		- Calcul de caractéristiques (intensité, volume, centroides).
///		- Filtrage morphologique basé sur les seuils d’intensité et de volume.
///		- Fitting gaussien (optionnel) pour affiner la position des objets.
///
/// Le résultat est un tableau de points, chacun contenant N_SEGMENTS doubles représentant
/// des informations telles que la position, l’intensité, les paramètres du fit, etc.
///
/// Le traitement est lancé via la méthode Processing().
/// </remarks>
class Palm
{
public:
	/// <summary>Constructeur de la classe <c>Palm</c>.</summary>
	///
	/// Initialise les structures internes pour effectuer une détection basée sur une décomposition en ondelettes et une segmentation intelligente.
	///
	/// <param name="image">Image d'entrée au format tableau 2D (grayscale, taille height × width).</param>
	/// <param name="points">Pointeur vers le tableau de sortie où seront stockés les points détectés (x, y, intensité, etc.).</param>
	/// <param name="pointsListSize">Taille totale allouée pour le tableau de points (doit être un multiple de N_SEGMENTS).</param>
	/// <param name="height">Hauteur de l'image (en pixels).</param>
	/// <param name="width">Largeur de l'image (en pixels).</param>
	/// <param name="thresholdVal">Seuil de détection appliqué sur les coefficients du plan d’ondelettes choisi.</param>
	/// <param name="watershedRatio">Ratio utilisé par l'algorithme de Watershed pour la séparation des objets.</param>
	/// <param name="gaussFit">Mode de fitting gaussien à appliquer :
	///		- 0 : aucun.
	///		- 1 : 2D avec sigma fixe.
	///		- 2 : 2D avec sigma variable.
	///		- 3 : 3D avec theta fixe.
	///		- 4 : 3D complet.</param>
	/// <param name="sigma">Valeur initiale de sigma utilisée dans le fitting (sigmaX = sigma et sigmaY = 2* sigma).</param>
	/// <param name="theta">Valeur initiale de l’orientation theta pour le fitting.</param>
	/// <param name="roi">Taille du ROI (fenêtre carrée centrée) utilisé pour le fitting autour de chaque point détecté.</param>
	Palm(const unsigned short* image, double* points, size_t pointsListSize, size_t height, size_t width,
	     double thresholdVal, double watershedRatio, size_t gaussFit, double sigma, double theta, size_t roi);

	/// <summary>Destructeur de la classe <c>Palm</c>.</summary>
	~Palm() = default;

	/// <summaryLance le traitement complet de détection sur l’image d’entrée.</summary>
	///
	/// Étapes principales :
	/// - Décomposition en ondelettes (ATrous).
	/// - Seuillage du plan choisi et création d'une carte binaire de segmentation.
	/// - Application d’un Watershed intelligent.
	/// - Calcul des caractéristiques des objets segmentés.
	/// - Application éventuelle de filtres morphologiques (volume, intensité).
	/// - Fitting gaussien optionnel autour des objets détectés.
	/// - Remplissage du tableau `points_` avec les résultats (coordonnées, intensité, sigma, etc.).
	///
	/// Le tableau `points_` est rempli avec des blocs de N_SEGMENTS double, représentant les paramètres pour chaque détection. Une ligne vide (valeurs = -1) est ajoutée à la fin.
	///
	/// <returns>Le nombre total de valeurs (doubles) écrites dans `points_`, y compris la ligne finale vide. (i.e. `nb_points * N_SEGMENTS + N_SEGMENTS`) </returns>
	int Processing();

private:
	const unsigned short* image_; ///< Image d'entrée (image 2D en niveaux de gris).
	double* points_;              ///< Tableau de points de sortie (rempli avec les résultats de détection).
	size_t height_,               ///< Hauteur de l'image.
	       width_,                ///< Largeur de l'image.
	       size_,                 ///< Taille totale de l'image (height × width).
	       nPointsMaxPerPlane_,   ///< Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).
	       currentPoint_,         ///< Indice du point actuel dans le tableau `points_`.
	       roi_;                  ///< Taille du ROI utilisé pour le fitting gaussien (fenêtre autour du point).
	double thresholdVal_,         ///< Seuil de détection utilisé après la décomposition en ondelettes.
	       watershedRatio_;       ///< Rapport utilisé dans l’algorithme de Watershed.

	size_t gaussFit_;   ///< Mode de fitting gaussien à appliquer (0 = aucun, 1 à 4 = variantes 2D ou 3D).
	double initSigmaX_, ///< Valeur initiale de sigmaX utilisée dans le fitting gaussien.
	       initSigmaY_, ///< Valeur initiale de sigmaY utilisée dans le fitting gaussien.
	       initTheta_;  ///< Valeur initiale de theta (orientation) utilisée dans le fitting gaussien.

	std::array<size_t, N_PLANES> operators_;                    ///< Opérateurs appliqués sur chaque plan d’ondelettes.

	ATrous aTrousPub_;                   ///< Objet interne pour la décomposition ondelettes (algorithme à trous).
	std::vector<unsigned short> segPub_; ///< Image segmentée (carte binaire des objets détectés).
	std::vector<double> outputPub_;      ///< Carte des valeurs utilisées pour la segmentation (après seuillage).
	Segmentation2D segmentPub_;          ///< Objet de segmentation intelligente 2D.
};
