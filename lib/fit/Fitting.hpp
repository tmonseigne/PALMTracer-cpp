///-------------------------------------------------------------------------------------------------
///
/// \file Fitting.hpp
/// \brief Classe mère pour l'ajustement.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <cstdint>
#include <vector>

#include "Constants.hpp"
#include "palm/Segmentation2D.hpp"

namespace palm_tracer {
/// <summary>Mode d'ajustement.</summary>
enum class fit_mode : uint8_t {
	none,              ///< Aucun Ajustement.
	gaussian_2d_fixed, ///< Ajustement Gaussien 2D avec Sigma fixé (uniquement X et Y).
	gaussian_2d_sigma, ///< Ajustement Gaussien 2D avec calcul du Sigma.
	gaussian_3d_fixed, ///< Ajustement Gaussien 3D avec Theta fixé et calcul du Sigma sur X et Y.
	gaussian_3d_theta, ///< Ajustement Gaussien 3D avec calcul du Sigma sur X et Y et du Theta.
	spline_3d_emccd,   ///< Ajustement par Spline pour Capteur EMCCD
	spline_3d_scmos    ///< Ajustement par Spline pour Capteur sCMOS
};

/// <summary>Classe d'ajustement gaussien.</summary>
class Fitting {
public:
	/// <summary>Constructeur de la classe <c>Fitting</c>. Initialise le nombre de paramètres du modèle en fonction du mode d'ajustement choisi.</summary>
	/// <param name="mode">Mode de fitting gaussien (2D, 2D avec sigma, 3D, 3D avec rotation).</param>
	/// <param name="image">Image 2D d'entrée (hauteur x largeur) en niveaux de gris.</param>
	/// <param name="height">Hauteur de l'image.</param>
	/// <param name="width">Largeur de l'image.</param>
	/// <param name="maxPoints">Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).</param>
	/// <param name="roi">Taille de la région d'intérêt autour du point.</param>
	/// <param name="sigmaX">Écart-type initial sur X.</param>
	/// <param name="sigmaY">Écart-type initial sur Y.</param>
	/// <param name="theta">Angle d'orientation (en radians) de la gaussienne si applicable.</param>
	explicit Fitting(fit_mode mode, const unsigned short* image, size_t height, size_t width, size_t maxPoints,
	                 size_t roi, double sigmaX, double sigmaY, double theta);

	/// <summary>Destructeur de la classe <c>Fitting</c>.</summary>
	virtual ~Fitting() = default;

	/// <summary>Exécute l'ajustement sur tous les points précédemments localisés par la segmentation.</summary>
	/// <param name="segmentation">Segmentation de l'image initiale. </param>
	/// <param name="points">Résultat de l'ajustement sur les points localisés.</param>
	void Fit(const Segmentation2D& segmentation, std::vector<std::array<double, N_LOCALIZATIONS>>& points);

protected:
	fit_mode mode_;                                    ///< Mode de fitting utilisé (ex. position seulement, ou position + orientation...).
	const unsigned short* image_;                      ///< Image 2D d'entrée (hauteur x largeur) en niveaux de gris.
	size_t n_,                                         ///< Nombre d'échantillons ou points de données à fitter.
	       nParams_,                                   ///< Nombre total de paramètres du modèle.
	       nParams2_,                                  ///< Nombre de paramètres au carré (utile pour les matrices jacobiennes ou hessiennes).
	       height_,                                    ///< Hauteur de l'image.
	       width_,                                     ///< Largeur de l'image.
	       maxPoints_,                                 ///< Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).
	       roi_,                                       ///< Taille de la région d'intérêt autour du point.
	       roiHalf_;                                   ///< Demi taille de la région d'intérêt autour du point.
	double sigmaX_,                                    ///< Écart-type initial sur X.
	       sigmaY_,                                    ///< Écart-type initial sur Y.
	       theta_;                                     ///< Angle d'orientation (en radians) de la gaussienne si applicable.
	static constexpr size_t N_ITER             = 50;   ///< Nombre maximal d'itérations autorisées dans l'algorithme d'optimisation.
	static constexpr double EPSILON_1          = 1e-6, ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2          = 1e-6, ///< Seuil de convergence pour le gradient.
	                        EPSILON_3          = 1e-5, ///< Seuil de convergence pour la variation des paramètres.
	                        EPSILON_4          = 1e-2; ///< Tolérance utilisée pour juger la stagnation ou instabilité.
	static constexpr double LAMBDA_0           = 1e-2, ///< Valeur initiale du paramètre lambda dans l'algorithme de Levenberg-Marquardt.
	                        LAMBDA_EPSILON_MAX = 1e7,  ///< Valeur maximale autorisée pour lambda.
	                        LAMBDA_EPSILON_MIN = 1e-7, ///< Valeur minimale autorisée pour lambda.
	                        LAMBDA_UP          = 11,   ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de lambda).
	                        LAMBDA_DN          = 9;    ///< Facteur de réduction de lambda lorsque la mise à jour est réussie.


	/// <summary>Exécute l'ajustement selon le mode sélectionné sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.
	/// L'algorithme utilise une méthode de Levenberg-Marquardt pour minimiser l'erreur entre le modèle et les données.</summary>
	/// <param name="infos">Informations sur le point localisé.</param>
	/// <param name="point">Résultat de l'ajustement sur le point localisé.</param>
	virtual void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) { fill(infos, point, infos.x_, infos.y_); }

	/// <summary></summary>
	/// <param name="infos">Informations initiales sur le point localisé.</param>
	/// <param name="point">Résultat de l'ajustement sur le point localisé.</param>
	/// <param name="x">Nouveau X.</param>
	/// <param name="y">Nouveau Y.</param>
	/// <param name="intensity">Nouvelle intensité.</param>
	/// <param name="sigmaX">Nouveau Sigma sur X.</param>
	/// <param name="sigmaY">Nouveau Sigma sur Y.</param>
	/// <param name="theta">Nouveau Theta.</param>
	/// <param name="rss">Nouveau RSS.</param>
	/// <param name="offset">Nouvel offset.</param>
	void fill(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point, double x, double y, double intensity = 0,
	          double sigmaX = 0, double sigmaY = 0, double theta = 0, double rss = 0, double offset = 0) const;
};
}
