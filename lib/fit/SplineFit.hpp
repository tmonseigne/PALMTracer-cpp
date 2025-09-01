///-------------------------------------------------------------------------------------------------
///
/// \file SplineFit.hpp
/// \brief Classe d'ajustement par spline.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <vector>

#include "Fitting.hpp"

namespace palm_tracer {
/// <summary>Classe d'ajustement par spline.</summary>
class SplineFit final : public Fitting {
public:
	explicit SplineFit(const fit_mode mode, const uint16_t* image, const Vec2I& dim, const size_t maxPoints,
	                   const double* fitParams) : Fitting(mode, image, dim, maxPoints, fitParams) { initialize(); }

protected:
	static constexpr double LAMBDA_0   = 1e-2, ///< Valeur initiale du paramètre lambda dans l'algorithme de Levenberg-Marquardt.
	                        LAMBDA_MIN = 1e-7, ///< Valeur minimale autorisée pour lambda.
	                        LAMBDA_MAX = 1e7,  ///< Valeur maximale autorisée pour lambda.
	                        LAMBDA_UP  = 11,   ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de lambda).
	                        LAMBDA_DN  = 9;    ///< Facteur de réduction de lambda lorsque la mise à jour est réussie.
	static constexpr double EPSILON_1  = 1e-6, ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2  = 1e-6, ///< Seuil de convergence pour le gradient.
	                        EPSILON_3  = 1e-5, ///< Seuil de convergence pour la variation des paramètres.
	                        EPSILON_4  = 1e-2; ///< Tolérance utilisée pour juger la stagnation ou instabilité.


	Vec3I splineDim_{};          ///< Dimensions de la spline (X, Y, Z)
	size_t splineSizeXY_{},      ///< Taille totale d’un plan XY. Utilisé pour accéder rapidement aux coefficients.
	       splineSize_{},        ///< Taille totale de l’ensemble du volume spline.
	       splineOffset_{};      ///< Offset utilisé pour indexer ou accéder efficacement aux sous-blocs ou aux patchs de la spline.
	double initZ_{},             ///< Position initiale en Z (origine physique ou de modélisation de la spline en profondeur).
	       deltaZ_{},            ///< Pas ou résolution entre deux couches Z successives dans la spline.
	       coeffRatio_{};        ///< Rapport d’intensité entre l’image et le modèle spline (pour corriger une mise à l’échelle ou un gain).
	std::vector<double> coeffs_; ///< Coefficients du modèle spline 3D, stockés de manière linéaire (ordre : [z][y][x]) pour un accès rapide.

	/// <summary>Parse les paramètres reçus pour l'ajustement.</summary>
	void parseParams() override;

	/// <summary>Exécute un ajustement de modèle spline sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const override;

	/// <summary>Calcule les dérivées du modèle spline 3D en un point donné (x, y, z).
	/// Remplit le vecteur des dérivées partielles (dudt) et évalue la valeur du modèle.</summary>
	/// <param name="pos">Position entière (X, Y, Z) dans la grille spline.</param>
	/// <param name="f">Valeurs de base spline pour l'intensité.</param>
	/// <param name="df">Référence vers les vecteurs contenant les dérivées de f par rapport à x (∂f/∂x), y (∂f/∂y) et z (∂f/∂z).</param>
	/// <param name="p">Vecteur des paramètres de la PSF (offset, intensité, position).</param>
	/// <param name="dudt">Vecteur des dérivées partielles du modèle par rapport aux paramètres.</param>
	/// <param name="model">>Valeur du modèle au point considéré.</param>
	void derivativeSpline(Vec3I& pos, const std::vector<double>& f, const std::array<std::vector<double>, 3>& df, const std::vector<double>& p,
	                      std::vector<double>& dudt, double& model) const;

	/// <summary>Calcule les contributions d'une base de fonctions splines 3D (valeurs et dérivées).
	/// Cette fonction calcule la valeur et les dérivées par rapport à x, y et z d'une fonction spline 3D basée sur des monômes (cx, cy, cz) à chaque point de l'espace (x, y, z).
	/// Le résultat est stocké dans trois tableaux (valeur, dérivées partielles) de taille SPLINE_PATCH_SIDE^3.
	/// </summary>
	/// <param name="x">Coordonnée x (locale, typiquement entre 0 et 1).</param>
	/// <param name="y">Coordonnée y (locale).</param>
	/// <param name="z">Coordonnée z (locale).</param>
	/// <param name="f">Référence vers le vecteur contenant les valeurs des fonctions de base.</param>
	/// <param name="df">Référence vers les vecteurs contenant les dérivées de f par rapport à x (∂f/∂x), y (∂f/∂y) et z (∂f/∂z).</param>
	void delta3D(const std::vector<double>& p, std::vector<double>& f, std::array<std::vector<double>, 3>& df, int& xi, int& yi, int& zi) const;

	/// <summary>Calcule la matrice alpha, le vecteur beta, et retourne le chi² (erreur) pour l'itération courante.
	/// Utilisé dans l'algorithme de Levenberg-Marquardt.</summary>
	/// <param name="grid">Vecteur des pixels de la zone d'intérêt : coordonnées (x, y) et valeur.</param>
	/// <param name="p">Paramètres du modèle</param>
	/// <returns>Valeur du χ², chi²  = r^T W r  (erreur quadratique totale)</returns>
	double lm(const std::vector<Pixel>& grid, const std::vector<double>& p, std::vector<double>& alpha, std::vector<double>& beta,
	          std::vector<double>& f, std::array<std::vector<double>, 3>& df, std::vector<double>& grad, int mode) const;


	[[nodiscard]] double gain(const std::vector<double>& deltaP, const std::vector<double>& alpha, const std::vector<double>& beta) const;
};
} // namespace palm_tracer
