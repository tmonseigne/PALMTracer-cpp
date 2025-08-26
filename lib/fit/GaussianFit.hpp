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
	static constexpr double LAMBDA_0   = 1e-2, ///< Valeur initiale du paramètre lambda dans l'algorithme de Levenberg-Marquardt.
	                        LAMBDA_MAX = 1e7,  ///< Valeur maximale autorisée pour lambda.
	                        LAMBDA_MIN = 1e-7, ///< Valeur minimale autorisée pour lambda.
	                        LAMBDA_UP  = 11,   ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de lambda).
	                        LAMBDA_DN  = 9;    ///< Facteur de réduction de lambda lorsque la mise à jour est réussie.
	static constexpr double EPSILON_1  = 1e-6, ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2  = 1e-6, ///< Seuil de convergence pour le gradient.
	                        EPSILON_3  = 1e-5, ///< Seuil de convergence pour la variation des paramètres.
	                        EPSILON_4  = 1e-2; ///< Tolérance utilisée pour juger la stagnation ou instabilité.

	double sigmaX_{1.0}, ///< Écart-type initial sur X.
	       sigmaY_{2.0}, ///< Écart-type initial sur Y.
	       theta_{0.0};  ///< Angle d'orientation (en radians) de la gaussienne si applicable.


	/// <summary>Parse les paramètres reçus pour l'ajustement.</summary>
	void parseParams() override;

	/// <summary>Exécute un ajustement de modèle Gaussien sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.
	/// L'algorithme utilise une méthode de Levenberg-Marquardt pour minimiser l'erreur entre le modèle et les données.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const override;

	//==================== Run Computes ====================
	/// <summary>Calcule le coefficient de corrélation au carré (R²) entre deux vecteurs normalisés.
	/// Cela donne une mesure de la qualité d'ajustement entre les données observées et le modèle.</summary>
	/// <param name="grid">Vecteur des pixels de la zone d'intérêt : coordonnées (x, y) et valeur (modifié en place).</param>
	/// <param name="hat">Vecteur des données ajustées (modifié en place).</param>
	/// <returns>Valeur de R², dans [0, 1].</returns>
	double computeR2(std::vector<Pixel>& grid, std::vector<double>& hat) const;

	/// <summary>Évalue la fonction PSF (Point Spread Function) en (x, y), selon le mode d'ajustement et les paramètres donnés.</summary>
	/// <param name="p">Vecteur des paramètres de la PSF (intensité, position, sigma, offset...).</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	/// <returns>Valeur de la PSF en (x, y).</returns>
	[[nodiscard]] double psf(const std::vector<double>& p, double x, double y) const;

	//==================== Matrix ====================
	/// <summary>Calcule le gradient de la PSF par rapport aux paramètres, pour un point (x, y), selon le mode d'ajustement.</summary>
	/// <param name="out">Pointeur vers le tableau de sortie du gradient (taille dépendant du mode).</param>
	/// <param name="p">Paramètres de la PSF.</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	void deriv(double* out, const std::vector<double>& p, double x, double y) const;

	/// <summary>Calcule la matrice alpha, le vecteur beta, et retourne le chi² (erreur) pour l'itération courante.
	/// Utilisé dans l'algorithme de Levenberg-Marquardt.</summary>
	/// <param name="grid">Vecteur des pixels de la zone d'intérêt : coordonnées (x, y) et valeur.</param>
	/// <param name="p">Paramètres du modèle</param>
	/// <param name="weight">Poids pour chaque point</param>
	/// <param name="hat">Vecteur des données ajustées (modifié en place), image "prédite" par le modèle au niveau de chaque pixel de la ROI.</param>
	/// <param name="alpha">Matrice de "courbure" Gauss-Newton, alpha = J^T W J</param>
	/// <param name="beta">Gradient du χ², beta  = J^T W r  (r = delta = résidus)</param>
	/// <returns>Valeur du χ², chi²  = r^T W r  (erreur quadratique totale)</returns>
	[[nodiscard]] double lm(const std::vector<Pixel>& grid, const std::vector<double>& p, const std::vector<double>& weight,
	                        std::vector<double>& hat, std::vector<double>& alpha, std::vector<double>& beta,
	                        std::vector<double>& delta, std::vector<double>& deltaW, std::vector<double>& grad) const;

	[[nodiscard]] double gain(const std::vector<double>& deltaP, const std::vector<double>& alpha, const std::vector<double>& beta) const;
};
} // namespace palm_tracer
