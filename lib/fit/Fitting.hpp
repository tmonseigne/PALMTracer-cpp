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
	spline_2d          ///< Ajustement par Spline
};

/// <summary>Classe d'ajustement gaussien.</summary>
class Fitting {
public:
	/// <summary>Constructeur de la classe <c>Fitting</c>. Initialise le nombre de paramètres du modèle en fonction du mode d'ajustement choisi.</summary>
	/// <param name="mode">Mode de fitting gaussien (2D, 2D avec sigma, 3D, 3D avec rotation).</param>
	/// <param name="image"></param>
	/// <param name="height"></param>
	/// <param name="width"></param>
	/// <param name="maxPoints"></param>
	explicit Fitting(fit_mode mode, const unsigned short* image, size_t height, size_t width, size_t maxPoints, size_t roi, double sigmaX, double sigmaY,
	                 double theta);

	/// <summary>Destructeur de la classe <c>Fitting</c>.</summary>
	virtual ~Fitting() = default;

	/// <summary</summary>
	/// <param name="segmentation"></param>
	/// <param name="points"></param>
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


	/// <summary>sdfs</summary>
	/// <param name="infos"></param>
	/// <param name="point"></param>
	virtual void fit(const SegmentInfo& infos, std::array<double, N_LOCALIZATIONS>& point) { fill(infos, point, infos.x_, infos.y_); }

	/// <summary></summary>
	/// <param name="infos"></param>
	/// <param name="point"></param>
	/// <param name="x"></param>
	/// <param name="y"></param>
	/// <param name="intensity"></param>
	/// <param name="sigmaX"></param>
	/// <param name="sigmaY"></param>
	/// <param name="theta"></param>
	/// <param name="rss"></param>
	/// <param name="offset"></param>
	void fill(const SegmentInfo& infos, std::array<double, N_LOCALIZATIONS>& point, double x, double y, double intensity = 0,
	          double sigmaX = 0, double sigmaY = 0, double theta = 0, double rss = 0, double offset = 0) const;
};
}
