///-------------------------------------------------------------------------------------------------
///
/// \file BSpline3D.hpp
/// \brief Représentation et traitement d'un volume 3D par régularisation B-spline de degré quelconque.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/07/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
/// \todo Attention le calcul final nécessite l'interpolation spline de matlab et il persiste une différence.
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <complex>
#include <vector>

#include "Structs.hpp"

namespace palm_tracer {
/// <summary>
/// Classe permettant de filtrer et interpoler un volume 3D selon un modèle B-spline régularisé.
/// Permet :
/// - le filtrage par pôles (Unser 1993) avec régularisation,
/// - le padding réfléchi sur chaque axe,
/// - l'interpolation régulière pour calcul de coefficients finaux d'interpolation spline continue.
/// </summary>
class BSpline3D {
public:
	static constexpr size_t N_DIMENSIONS = 3; ///< Nombre de dimensions (X, Y, Z).

	Coord3D dim_{},       ///< Dimensions d'origine du volume (sans padding).
	        dimPadded_{}, ///< Dimensions après padding réfléchi.
	        dimFinal_{};  ///< Dimensions des coefficients finaux (dim - 1).

	size_t degree_ = 3;                         ///< Degré de la B-spline (entre 0 et 7).
	std::array<double, N_DIMENSIONS> lambda_{}; ///< Facteur de régularisation (0 = pas de régularisation, ≠ 0 autorisé uniquement pour les degrés impairs).

	std::array<std::vector<std::complex<double>>, N_DIMENSIONS> filters_{}; ///< Liste des pôles par dimension. (réels ou complexes conjugués)
	std::array<double, N_DIMENSIONS> c0_{};                                 ///< Coefficient multiplicatif final à appliquer au signal par dimension.

	std::vector<double> coeffs_,       ///< Coefficients après filtrage.
	                    coeffsPadded_, ///< Coefficients avec padding réfléchi.
	                    coeffsHDR_,    ///< Volume interpolé à haute résolution (matlab-like).
	                    coeffsFinal_;  ///< Coefficients finaux pour interpolation continue.

	/// <summary>Constructeur par défaut de la classe <c>BSpline3D</c>.</summary>
	BSpline3D() = default;

	/// <summary>Constructeur de la classe <c>BSpline3D</c>: calcule et filtre les coefficients B-spline d'un volume 3D.</summary>
	/// <param name="volume">Volume 3D d'entrée (vectorisé en 1D).</param>
	/// <param name="dim">Dimensions du volume (X, Y, Z).</param>
	/// <param name="degree">Degré de la B-spline (0 à 7).</param>
	/// <param name="lambda">Paramètres de régularisation pour chaque dimension.</param>
	BSpline3D(const std::vector<double>& volume, const Coord3D& dim, size_t degree = 3, const Vec3& lambda = {});

	/// <summary>Destructeur de la classe <c>BSpline3D</c>.</summary>
	~BSpline3D() = default;

	/// <summary> Évalue la fonction de base B-spline de degré donné sur un vecteur de positions.</summary>
	/// <param name="x">position (peut contenir des négatifs, réels, etc).</param>
	/// <param name="deg">Degré de la B-spline (0 à 7 pris en charge explicitement).</param>
	/// <returns>Valeur de la base B-spline à la position x.</returns>
	static double Eval(double x, size_t deg);

protected:
	/// <summary>Récupère les coefficients finaux d'interpolation continue (64 par voxel).
	/// Nécessite l'interpolation triplée du volume et une inversion matricielle.</summary>
	void getFinalCoeff();

	/// <summary>Initialise les filtres B-spline selon le degré et la régularisation choisie.</summary>
	void getPolesFilter();

	/// <summary>Applique un filtre exponentiel symétrique (Unser, 1993) sur l'axe spécifié.</summary>
	/// <param name="axis">Axe concerné (0=X, 1=Y, 2=Z).</param>
	/// <param name="poles">Pôles à appliquer.</param>
	/// <param name="gain">Gain global à appliquer (fonction de lambda).</param>
	void applyFilter(size_t axis, const std::vector<std::complex<double>>& poles, double gain);

	/// <summary>Interpolant régulier du volume régularisé (utilise `interpolate_at` sur tous les voxels).</summary>
	void interpolateHDR();

	/// <summary>Interpolation à une position réelle dans le volume paddé (spline continue).</summary>
	/// <param name="x">Coordonnée X (réelle).</param>
	/// <param name="y">Coordonnée Y (réelle).</param>
	/// <param name="z">Coordonnée Z (réelle).</param>
	/// <returns>Valeur interpolée au point (x, y, z).</returns>
	[[nodiscard]] double interpolate_at(double x, double y, double z) const;

	/// <summary>Résolution d'un polynôme réel (3e, 5e ou 7e degré) par Cardano ou Durand-Kerner.</summary>
	/// <param name="coeffs">Coefficients du polynôme (en partant du degré le plus haut).</param>
	/// <param name="tol">Tolérance pour la convergence (Durand-Kerner).</param>
	/// <param name="max_iter">Nombre maximal d'itérations.</param>
	/// <returns>Vecteur des racines complexes.</returns>
	[[nodiscard]] static std::vector<std::complex<double>> solvePolynomial(const std::vector<double>& coeffs, double tol = 1e-12, int max_iter = 100);

	/// <summary> Applique un filtre exponentiel symétrique (causal/anticausal) (Unser, 1993) à un vecteur 1D.</summary>
	/// <param name="x">Entrée (vecteur complexe).</param>
	/// <param name="c0">Constante de normalisation (dépend du pôle).</param>
	/// <param name="zi">Pôle du filtre (peut être complexe).</param>
	/// <param name="k0">Profondeur de réflexion initiale.</param>
	/// <param name="kVec">Vecteur d'indices réfléchis pour initialisation (taille K0).</param>
	/// <returns>Signal filtré (même taille que `x`).</returns>
	static std::vector<std::complex<double>> symExpFilt(const std::vector<std::complex<double>>& x, const std::complex<double>& c0,
	                                                    const std::complex<double>& zi, size_t k0, const std::vector<size_t>& kVec);
};
} // namespace palm_tracer
