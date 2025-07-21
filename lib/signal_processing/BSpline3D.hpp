///-------------------------------------------------------------------------------------------------
///
/// \file Spline.hpp
/// \brief Classe de calibration de la spline pour l'ajustement.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/07/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <complex>
#include <vector>

#include "Structs.hpp"

namespace palm_tracer {
/// <summary>Classe de calibration de la spline pour l'ajustement.</summary>
class BSpline3D {
public:
	static constexpr size_t N_DIMENSIONS = 3;

	Coord3D dim_{}, dimCoeffs_{}, dimFinal_{};

	size_t degree_;                             ///< Degré de la B-spline (0 à 7)
	std::array<double, N_DIMENSIONS> lambda_{}; ///< Facteur de régularisation (≠ 0 autorisé uniquement pour les degrés impairs).

	std::array<std::vector<std::complex<double>>, N_DIMENSIONS> filter_{}; ///< Vecteur de pôles du filtre (réels ou complexes conjugués)
	std::array<double, N_DIMENSIONS> c0_{};                                ///< Coefficient multiplicatif final à appliquer au signal.

	std::vector<double> coeffs_, coeffsHDR_, coeffsFinal_;

	BSpline3D(const std::vector<double>& volume, const Coord3D& dim, size_t degree = 3, const Vec3& lambda = {});
	~BSpline3D() = default;

protected:
	void getFinalCoeff();

	/// <summary> Calcule les pôles du filtre B-spline pour une régularisation directe.</summary>
	void getPolesFilter();

	/// <summary>Applies the filter.</summary>
	/// <param name="axis">The axis.</param>
	/// <param name="poles">The poles.</param>
	/// <param name="gain">The gain.</param>
	void applyFilter(size_t axis, const std::vector<std::complex<double>>& poles, double gain);

	void interpolateHDR();

	[[nodiscard]] double interpolate_at(double x, double y, double z) const;

	/// <summary>Résout un polynôme réel à coefficients double de degré 3, 5 ou 7. Retourne les racines complexes.</summary>
	/// <param name="coeffs">The coeffs.</param>
	/// <param name="tol">The tol.</param>
	/// <param name="max_iter">The maximum iter.</param>
	/// <returns></returns>
	[[nodiscard]] static std::vector<std::complex<double>> solvePolynomial(const std::vector<double>& coeffs, double tol = 1e-12, int max_iter = 100);

	/// <summary> Applique un filtre exponentiel symétrique (Unser, 1993) à un vecteur 1D. </summary>
	/// <param name="x">Entrée brute (non modifiée).</param>
	/// <param name="c0">Constante de symétrisation (dépend du pôle).</param>
	/// <param name="zi">Pôle du filtre (peut être complexe).</param>
	/// <param name="k0">Profondeur initiale de réflexion.</param>
	/// <param name="kVec">Indices réfléchis pour initialisation (taille K0).</param>
	/// <returns>Vecteur filtré (même taille que X).</returns>
	static std::vector<std::complex<double>> symExpFilt(const std::vector<std::complex<double>>& x, const std::complex<double>& c0,
	                                                    const std::complex<double>& zi, size_t k0, const std::vector<size_t>& kVec);

	/// <summary> Évalue la fonction de base B-spline de degré donné sur un vecteur de positions. </summary>
	/// <param name="x">position (peut contenir des négatifs, réels, etc).</param>
	/// <param name="deg">Degré de la B-spline (0 à 7 pris en charge explicitement).</param>
	/// <returns>Valeur de la base B-spline évaluée.</returns>
	static double eval(double x, size_t deg);
};
} // namespace palm_tracer
