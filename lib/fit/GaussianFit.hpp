///-------------------------------------------------------------------------------------------------
///
/// \file GaussianFit.hpp
/// \brief Classe d'ajustement gaussien.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <vector>

#include "Fitting.hpp"

namespace palm_tracer {
/// <summary>Classe d'ajustement gaussien.</summary>
class GaussianFit final : public Fitting {
public:
	explicit GaussianFit(const fit_mode mode, const uint16_t* image, const Vec2I& dim, const size_t maxPoints,
	                     const double* fitParams) : Fitting(mode, image, dim, maxPoints, fitParams) { initialize(); }

protected:
	static constexpr double LAMBDA_0   = 1e-2, ///< Valeur initiale du paramètre λ dans l'algorithme de Levenberg-Marquardt.
	                        LAMBDA_MIN = 1e-7, ///< Valeur minimale autorisée pour λ.
	                        LAMBDA_MAX = 1e7,  ///< Valeur maximale autorisée pour λ.
	                        LAMBDA_UP  = 11,   ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de λ).
	                        LAMBDA_DN  = 9;    ///< Facteur de réduction de λ lorsque la mise à jour est réussie.
	static constexpr double EPSILON_1  = 1e-6, ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2  = 1e-6, ///< Seuil de convergence pour la variation des paramètres.
	                        EPSILON_3  = 1e-5, ///< Seuil de convergence pour l'erreur moyenne par pixel.
	                        EPSILON_4  = 1e-2; ///< Tolérance utilisée pour juger la stagnation ou instabilité.

	double sigmaX_{1.0}, ///< Écart-type initial sur X.
	       sigmaY_{2.0}, ///< Écart-type initial sur Y.
	       theta_{0.0};  ///< Angle d'orientation (en radians) de la gaussienne si applicable.


	/// <summary>Ensemble des objets nécessaires au déroulement de l’algorithme de Levenberg-Marquardt.
	/// Chaque membre correspond à une variable intermédiaire utilisée dans les calculs d’ajustement non-linéaire.</summary>
	struct LM_Workspace {
		std::vector<double> alpha,            ///< Matrice de courbure (Gauss-Newton), de taille np² : α = Jᵀ × W × J.
		                    alphaR,           ///< Copie d'alpha avec régularisation diagonale par λ (αR = α + λ diag(α)).
		                    alphaRi,          ///< Inverse d'alphaR  (obtenue par inversion de Gauss-Jordan).
		                    beta,             ///< Gradient du χ² (β = Jᵀ × W × r), taille np.
		                    delta,            ///< Résidus non pondérés : différence entre données observées et modèle (δ = y - ŷ).
		                    deltaW,           ///< Résidus pondérés : δW = W × δ.
		                    deltaP,           ///< Incrément de mise à jour des paramètres (δp), calculé à chaque itération.
		                    lBetaDeltaP,      ///< Terme auxiliaire : λ × δp + β, utilisé pour calculer le gain prédictif.
		                    grad,             ///< Gradient de la PSF par rapport aux paramètres, évalué pixel par pixel.
		                    pTry,             ///< Paramètres candidats (p + δp) après mise à jour et clampés aux bornes.
		                    hat,              ///< Valeurs prédites par le modèle (ŷ) pour chaque pixel de la ROI.
		                    weight;           ///< Poids associés à chaque pixel (permet d’ajuster l’influence des points).
		std::vector<Pixel> grid;              ///< Liste des pixels de la ROI, contenant coordonnées (x, y) et intensité mesurée.
		std::array<std::vector<double>, 3> p; ///< Paramètres du modèle : borne inférieure, valeur courante et borne supérieure de chaque paramètre.

		/// <summary>Redimensionne et initialise les différents vecteurs nécessaires à l’algorithme de Levenberg-Marquardt.</summary>
		/// <param name="n">Nombre de pixels de la ROI (taille de grid, hat, delta...).</param>
		/// <param name="np">Nombre de paramètres du modèle.</param>
		/// <param name="np2">Nombre de paramètres au carré (np²).</param>
		void Resize(const size_t n, const size_t np, const size_t np2)
		{
			alpha.resize(np2);
			alphaR.resize(np2);
			alphaRi.resize(np2);

			beta.resize(np);
			grad.resize(np);
			deltaP.resize(np);
			lBetaDeltaP.resize(np);
			pTry.resize(np);

			hat.resize(n);
			weight.resize(n, 1);
			delta.resize(n);
			deltaW.resize(n);
			grid.resize(n);
		}
	};

	//==================== Common Functions ====================
	/// <summary>Parse les paramètres reçus pour l'ajustement.</summary>
	void parseParams() override;

	/// <summary>Exécute un ajustement de modèle Gaussien sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.
	/// L'algorithme utilise une méthode de Levenberg-Marquardt pour minimiser l'erreur entre le modèle et les données.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const override;

	//==================== Levenberg-Marquardt ====================
	/// <summary>Calcule la ROI du modèle, la matrice alpha, le vecteur beta, et retourne l'erreur (chi²) pour l'itération courante.
	/// Utilisé dans l'algorithme de Levenberg-Marquardt.</summary>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme de Levenberg-Marquardt.</param>
	/// <returns>Valeur du χ², chi² = r^T × W × r  (erreur quadratique totale)</returns>
	[[nodiscard]] double lm(LM_Workspace& ws) const;

	/// <summary>Estime le gain de l'iteration acutelle (2 × δPᵀ × β - δPᵀ × α × δP).</summary>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme de Levenberg-Marquardt.</param>
	/// <returns>Gain estimé.</returns>
	[[nodiscard]] double gain(const LM_Workspace& ws) const;

	/// <summary>Calcul l'erreur (chi²) entre données observées et modèle candidat.</summary>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme de Levenberg-Marquardt.</param>
	/// <returns>Erreur calculée.</returns>
	[[nodiscard]] double getError(const LM_Workspace& ws) const;

	/// <summary>Calcule le coefficient de corrélation au carré (R²) entre données observées et modèle.
	/// Cela donne une mesure de la qualité d'ajustement entre les données observées et le modèle.</summary>
	/// <param name="ws">Ensemble des objets nécessaires au déroulement de l’algorithme de Levenberg-Marquardt.</param>
	/// <returns>Valeur de R², dans [0, 1].</returns>
	[[nodiscard]] double computeR2(LM_Workspace& ws) const;

	//==================== Run Computes ====================
	/// <summary>Évalue la fonction PSF (Point Spread Function) en (x, y), selon le mode d'ajustement et les paramètres donnés.</summary>
	/// <param name="p">Vecteur des paramètres de la PSF (intensité, position, sigma, offset...).</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	/// <returns>Valeur de la PSF en (x, y).</returns>
	[[nodiscard]] double psf(const std::vector<double>& p, double x, double y) const;

	/// <summary>Calcule le gradient de la PSF par rapport aux paramètres, pour un point (x, y), selon le mode d'ajustement.</summary>
	/// <param name="grad">Vecteur de sortie du gradient (taille dépendant du mode).</param>
	/// <param name="p">Paramètres de la PSF.</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	void deriv(std::vector<double>& grad, const std::vector<double>& p, double x, double y) const;
};
} // namespace palm_tracer
