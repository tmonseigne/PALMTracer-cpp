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
#include <array>
#include <vector>

#include "Constants.hpp"
#include "Fitting.hpp"

namespace palm_tracer {
/// <summary>Classe d'ajustement par spline.</summary>
class SplineFit final : public Fitting {
public:
	explicit SplineFit(const fit_mode mode, const uint16_t* image, const Vec2I& dim, const size_t maxPoints,
	                   const double* fitParams) : Fitting(mode, image, dim, maxPoints, fitParams) { initialize(); }

protected:
	static constexpr double LAMBDA_0   = 0.1,   ///< Valeur initiale du paramètre λ dans l'algorithme du Maximum de vraisemblance.
	                        LAMBDA_MIN = 1e-7,  ///< Valeur minimale autorisée pour λ.
	                        LAMBDA_MAX = 1e7,   ///< Valeur maximale autorisée pour λ.
	                        LAMBDA_UP  = 11,    ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de λ).
	                        LAMBDA_DN  = 9;     ///< Facteur de réduction de λ lorsque la mise à jour est réussie.
	static constexpr double EPSILON_1  = 1e-6,  ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2  = 1e-10; ///< avant calcul de la matrice inverse

	Vec3I splineDim_{};          ///< Dimensions de la spline (X, Y, Z)
	size_t splineSizeXY_{},      ///< Taille totale d’un plan XY. Utilisé pour accéder rapidement aux coefficients.
	       splineSize_{},        ///< Taille totale de l’ensemble du volume spline.
	       splineOffset_{};      ///< Offset utilisé pour indexer ou accéder efficacement aux sous-blocs ou aux patchs de la spline.
	double initZ_{},             ///< Position initiale en Z (origine physique ou de modélisation de la spline en profondeur).
	       deltaZ_{},            ///< Pas ou résolution entre deux couches Z successives dans la spline.
	       coeffRatio_{};    ///< Rapport d’intensité entre l’image et le modèle spline (pour corriger une mise à l’échelle ou un gain).
	std::vector<double> coeffs_; ///< Coefficients du modèle spline 3D, stockés de manière linéaire (ordre : [z][y][x]) pour un accès rapide.

	/// <summary>Ensemble des objets nécessaires au déroulement de l’algorithme du Maximum de vraisemblance (Maximum likelihood estimation).
	/// Chaque membre correspond à une variable intermédiaire utilisée dans les calculs d’ajustement non-linéaire.</summary>
	struct MLE_Workspace {
		std::vector<double> jacobian,                       ///<
		                    hessian,                        ///<
		                    l,                              ///< Matrice triangulaire inférieure (L) du système linéaire L × U × x = b.
		                    u,                              ///< Matrice triangulaire supérieure (U) du système linéaire L × U × x = b.
		                    dudt,                           ///< Vecteur des dérivées partielles du modèle par rapport aux paramètres.
		                    jump,                           ///<
		                    update,                         ///<
		                    p_try,                          ///< Paramètres candidats (p + δp formule différente ici) après mise à jour et clampés aux bornes.
		                    f;                              ///<
		std::vector<Pixel> grid;                            ///<
		std::array<std::vector<double>, 3> p{};             ///< Valeurs de base spline pour l'intensité.
		std::array<std::vector<double>, 3> df{};            ///< Dérivées de f par rapport à x (∂f/∂x), y (∂f/∂y) et z (∂f/∂z).
		double lambda = LAMBDA_0,                           ///<
		       model  = 0,                                  ///< Valeur du modèle au point considéré.
		       err    = std::numeric_limits<double>::max(), ///<
		       mu     = 1 + lambda;                         ///<


		/// <summary>Redimensionne et initialise les différents vecteurs nécessaires à l’algorithme du Maximum de vraisemblance.</summary>
		/// <param name="n">Nombre de pixels de la ROI.</param>
		/// <param name="np">Nombre de paramètres du modèle.</param>
		/// <param name="np2">Nombre de paramètres au carré (np²).</param>
		void Resize(const size_t n, const size_t np, const size_t np2)
		{
			hessian.resize(np2);
			l.resize(np2);
			u.resize(np2);

			jacobian.resize(np);
			dudt.resize(np);
			jump.resize(np);
			update.resize(np);
			p_try.resize(np);

			grid.resize(n);

			f.resize(SPLINE_PATCH_SIZE);
			for (auto& d : df) { d.resize(SPLINE_PATCH_SIZE); }
		}
	};

	//==================== Common Functions ====================
	/// <summary>Parse les paramètres reçus pour l'ajustement.</summary>
	void parseParams() override;

	/// <summary>Exécute un ajustement de modèle spline sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const override;

	//==================== Maximum de vraisemblance ====================
	/// <summary>Met à jour le jacobien et l'hessienne pour le modèle MLE (Maximum Likelihood Estimation).
	/// Utilise une vraisemblance basée sur la statistique de Poisson.</summary>
	/// <param name="data">Valeur observée dans l'image.</param>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme du Maximum de vraisemblance.</param>
	void updateJacobianHessian(double data, MLE_Workspace& ws) const;

	/// <summary>Met à jour les paramètres d’ajustement de type Levenberg-Marquardt.
	/// Modifie mu et lambda selon si l’itération précédente était un succès.</summary>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme du Maximum de vraisemblance.</param>
	/// <param name="up">True si l’itération précédente a échoué (augmentation du lambda).</param>
	static void updateMuLambda(MLE_Workspace& ws, bool up);

	/// <summary>Résout un système linéaire LU x = b avec matrices l (lower) et u (upper).
	/// Effectue une substitution avant (Ly = b) puis arrière (Ux = y).</summary>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme du Maximum de vraisemblance.</param>
	void luEvaluate(MLE_Workspace& ws) const;

	//==================== Run Computes ====================
	/// <summary>Calcule les dérivées du modèle spline 3D en un point donné (x, y, z).
	/// Remplit le vecteur des dérivées partielles (dudt) et évalue la valeur du modèle.</summary>
	/// <param name="pos">Position entière (X, Y, Z) dans la grille spline.</param>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme du Maximum de vraisemblance.</param>
	void derivativeSpline(Vec3I& pos, MLE_Workspace& ws) const;

	/// <summary>Calcule les contributions d'une base de fonctions splines 3D (valeurs et dérivées).
	/// Cette fonction calcule la valeur et les dérivées par rapport à x, y et z d'une fonction spline 3D basée sur des monômes (cx, cy, cz) à chaque point de l'espace (x, y, z).
	/// Le résultat est stocké dans trois tableaux (valeur, dérivées partielles) de taille SPLINE_PATCH_SIDE^3.
	/// </summary>
	/// <param name="pos">Coordonnées (x, y, z) (locale, typiquement entre 0 et 1).</param>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme du Maximum de vraisemblance.</param>
	static void delta3D(const Vec3D& pos, MLE_Workspace& ws);

	/// <summary>Calcul un coefficient multiplicateur sur la base du centre du voxel central du volume.</summary>
	double computeCoeffRatio() const;
};
} // namespace palm_tracer
