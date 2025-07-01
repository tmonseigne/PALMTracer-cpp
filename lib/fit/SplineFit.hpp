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
	explicit SplineFit(const fit_mode mode, const uint16_t* image, const size_t height, const size_t width, const size_t maxPoints,
	                   const double* fitParams) : Fitting(mode, image, height, width, maxPoints, fitParams) { initialize(); }

protected:
	static constexpr double LAMBDA_0   = 0.1,   ///< Valeur initiale du paramètre lambda dans l'algorithme de Levenberg-Marquardt.
	                        LAMBDA_MAX = 1e7,   ///< Valeur maximale autorisée pour lambda.
	                        LAMBDA_MIN = 1e-7,  ///< Valeur minimale autorisée pour lambda.
	                        LAMBDA_UP  = 10,    ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de lambda).
	                        LAMBDA_DN  = 0.1;   ///< Facteur de réduction de lambda lorsque la mise à jour est réussie.
	static constexpr double EPSILON_1  = 1e-6,  ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2  = 1e-10; ///< avant calcul de la matrice inverse

	size_t splineSizeX_{},       ///< Taille de la spline en X (nombre de points dans la direction horizontale).
	       splineSizeY_{},       ///< Taille de la spline en Y (nombre de points dans la direction verticale).
	       splineSizeZ_{},       ///< Taille de la spline en Z (nombre de plans en profondeur).
	       splineSizeXY_{},      ///< Taille totale d’un plan XY (splineSizeX_ * splineSizeY_). Utilisé pour accéder rapidement aux coefficients.
	       splineSizeXYZ_{},     ///< Taille totale de l’ensemble du volume spline (splineSizeX_ * splineSizeY_ * splineSizeZ_).
	       splineOffset_{};      ///< Offset utilisé pour indexer ou accéder efficacement aux sous-blocs ou aux patchs de la spline.
	double initZ_{},             ///< Position initiale en Z (origine physique ou de modélisation de la spline en profondeur).
	       deltaZ_{},            ///< Pas ou résolution entre deux couches Z successives dans la spline.
	       intensityRatio_{};    ///< Rapport d’intensité entre l’image et le modèle spline (pour corriger une mise à l’échelle ou un gain).
	std::vector<double> coeffs_; ///< Coefficients du modèle spline 3D, stockés de manière linéaire (ordre : [z][y][x]) pour un accès rapide.

	/// <summary>Parse les paramètres reçus pour l'ajustement.</summary>
	void parseParams() override;

	/// <summary>Exécute un ajustement de modèle spline sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) override;

	/// <summary>Calcule les dérivées du modèle spline 3D en un point donné (x, y, z).
	/// Remplit le vecteur des dérivées partielles (dudt) et évalue la valeur du modèle.</summary>
	/// <param name="x">Indice x dans la grille spline (colonne).</param>
	/// <param name="y">Indice y dans la grille spline (ligne).</param>
	/// <param name="z">Indice z dans la grille spline (profondeur).</param>
	/// <param name="f">Valeurs de base spline pour l'intensité.</param>
	/// <param name="dxf">Dérivées de f par rapport à x (∂f/∂x).</param>
	/// <param name="dyf">Dérivées de f par rapport à y (∂f/∂y).</param>
	/// <param name="dzf">Dérivées de f par rapport à z (∂f/∂z).</param>
	/// <param name="p">Vecteur des paramètres de la PSF (offset, intensité, position).</param>
	/// <param name="dudt">Vecteur des dérivées partielles du modèle par rapport aux paramètres.</param>
	/// <param name="model">>Valeur du modèle au point considéré.</param>
	void derivativeSpline(size_t x, size_t y, size_t z, const std::vector<double>& f, const std::vector<double>& dxf, const std::vector<double>& dyf,
	                      const std::vector<double>& dzf, const std::vector<double>& p, std::vector<double>& dudt, double& model) const;

	/// <summary> Met à jour le jacobien et l'hessienne pour le modèle MLE (Maximum Likelihood Estimation).
	/// Utilise une vraisemblance basée sur la statistique de Poisson. </summary>
	/// <param name="data">Valeur observée dans l'image.</param>
	/// <param name="model">Valeur prédite par le modèle.</param>
	/// <param name="err">Erreur totale accumulée (log-vraisemblance).</param>
	/// <param name="jacobian">Vecteur du gradient (jacobien).</param>
	/// <param name="hessian">Matrice Hessienne (seconde dérivée).</param>
	/// <param name="dudt">Vecteur des dérivées partielles du modèle.</param>
	void updateJacobianHessian(double data, double model, double& err, std::vector<double>& jacobian,
	                           std::vector<double>& hessian, const std::vector<double>& dudt) const;

	/// <summary> Met à jour les paramètres d’ajustement de type Levenberg-Marquardt.
	/// Modifie mu et lambda selon si l’itération précédente était un succès.</summary>
	/// <param name="mu">Facteur de régularisation.</param>
	/// <param name="lambda">Paramètre de contrôle LM.</param>
	/// <param name="up">True si l’itération précédente a échoué (augmentation du lambda).</param>
	static void updateMuLambda(double& mu, double& lambda, bool up);

	/// <summary>Résout un système linéaire LU x = b avec matrices l (lower) et u (upper).
	/// Effectue une substitution avant (Ly = b) puis arrière (Ux = y).</summary>
	/// <param name="l">Matrice triangulaire inférieure (L).</param>
	/// <param name="u">Matrice triangulaire supérieure (U).</param>
	/// <param name="b">Vecteur second membre (b).</param>
	/// <param name="x">Vecteur solution du système (x).</param>
	void luEvaluate(const std::vector<double>& l, const std::vector<double>& u, const std::vector<double>& b, std::vector<double>& x) const;

	/// <summary>Calcule les contributions d'une base de fonctions splines 3D (valeurs et dérivées).
	/// Cette fonction calcule la valeur et les dérivées par rapport à x, y et z d'une fonction spline 3D basée sur des monômes (cx, cy, cz) à chaque point de l'espace (x, y, z).
	/// Le résultat est stocké dans trois tableaux (valeur, dérivées partielles) de taille SPLINE_PATCH_SIDE^3.
	/// </summary>
	/// <param name="x">Coordonnée x (locale, typiquement entre 0 et 1).</param>
	/// <param name="y">Coordonnée y (locale).</param>
	/// <param name="z">Coordonnée z (locale).</param>
	/// <param name="f">Référence vers le vecteur contenant les valeurs des fonctions de base.</param>
	/// <param name="dxf">Référence vers le vecteur contenant les dérivées de f par rapport à x (∂f/∂x).</param>
	/// <param name="dyf">Référence vers le vecteur contenant les dérivées de f par rapport à y (∂f/∂y).</param>
	/// <param name="dzf">Référence vers le vecteur contenant les dérivées de f par rapport à z (∂f/∂z).</param>
	static void delta3D(double x, double y, double z, std::vector<double>& f, std::vector<double>& dxf, std::vector<double>& dyf, std::vector<double>& dzf);
};
} // namespace palm_tracer
