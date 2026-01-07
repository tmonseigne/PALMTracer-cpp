///-------------------------------------------------------------------------------------------------
///
/// \file AstigFit.hpp
/// \brief Classe d'ajustement non linéaire (Levenberg–Marquardt) d'une courbe d'astigmatisme \f$\sigma(z)\f$.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 05/01/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <vector>

#include "Structs.hpp"

namespace palm_tracer::fit {
/// <summary>Coefficients d'astigmatisme complets [0]=X, [1]=Y. Les coefficients sont : [z0,w,c3,c4,A].
///
/// - **z0** : Décalage axial (z0) en nm.
/// - **w** : Échelle axiale (w) en nm. Doit être strictement supérieur à 0.
/// - **c3** : Terme cubique sans dimension.
/// - **c4** : Terme quartique sans dimension.
/// - **A** : Amplitude (sigma au foyer) en nm.
/// </summary>
using astig_model = std::array<std::array<double, 5>, 2>; // [0]=X,[1]=Y ; [z0,w,c3,c4,A]

/// <summary>Ajustement non linéaire (Levenberg–Marquardt) d'une courbe d'astigmatisme \f$\sigma(z)\f$, et estimation de Z à partir d'une paire \f$(\sigma_x, \sigma_y) \f$.</summary>
/// <remarks>
/// Hypothèses importantes :
/// - Les entrées \f$\sigma\f$ sont en pixels.
/// - Les entrées z sont en nanomètres.
/// - Les paramètres retournés <c>A</c> sont en nanomètres.
/// - <see cref="sampling_"/> et <see cref="pixelSize_"/> pilotent la conversion nm ↔ px.
/// </remarks>
class AstigFit {
public:
	/// <summary>Constructeur de la classe <see cref="AstigFit"/>.</summary>
	/// <param name="pixelSize">Taille de pixel (nm/px). Doit être strictement > 0.</param>
	/// <param name="sampling">Facteur de sur-échantillonnage (sans dimension). 1 = pixels natifs caméra ; 2 = reconstruction sur grille 2× plus fine, etc.</param>
	/// <remarks>Si le pipeline n'upsample pas les coordonnées / sigmas, laisser <paramref name="sampling"/> à 1.</remarks>
	explicit AstigFit(const double pixelSize = 160, const double sampling = 1.0) : pixelSize_(pixelSize), sampling_(sampling) { }

	/// <summary>Calibre (fit LM) les deux courbes \f$\sigma_x(z)\f$ et \f$\sigma_y(z)\f$ à partir d'un tableau (σx, σy, z).</summary>
	/// <param name="input">
	/// Points de calibration, taille N :
	/// - <c>input[i].x</c> = \f$\sigma_x\f$ en pixels
	/// - <c>input[i].y</c> = \f$\sigma_y\f$ en pixels
	/// - <c>input[i].z</c> = z en nanomètres
	/// </param>
	/// <returns>
	/// Le modèle 2×5 :
	/// - <c>res[0] = [z0x, wx, c3x, c4x, Ax]</c>
	/// - <c>res[1] = [z0y, wy, c3y, c4y, Ay]</c>
	/// </returns>
	/// <remarks>
	/// - Les axes X et Y sont ajustés indépendamment (<see cref="fitAxis"/>).
	/// - En cas d'entrée invalide (tailles incohérentes, N trop petit), l'axe concerné renvoie un vecteur rempli de -1.
	/// </remarks>
	[[nodiscard]] astig_model Calibrate(const std::vector<Vec3D>& input) const;

	/// <summary>Estime z (nm) en minimisant l'erreur quadratique entre \f$(\sigma_x, \sigma_y)\f$ mesurés et prédits par le modèle.</summary>
	/// <param name="model">Modèle d'astigmatisme (résultat de <see cref="Calibrate"/>).</param>
	/// <param name="sigma">Mesure \f$(\sigma_x, \sigma_y)\f$ en pixels.</param>
	/// <param name="zLimit">Intervalle de recherche initial : <c>[zMin, zMax]</c> en nm.</param>
	/// <param name="zStep">Nombre d'échantillons pour le pré-scan (grille grossière) avant raffinement.
	/// Doit être ≥ 3 (sinon division par <c>zStep-1</c>) et impaire si l'on veut le point 0 dans les choix initiaux.</param>
	/// <returns>Z (nm) minimisant l'erreur sur l'intervalle.</returns>
	/// <remarks>
	/// Stratégie :
	/// 1) pré-scan sur <paramref name="zLimit"/> pour trouver un bon point initial,
	/// 2) raffinement par <see cref="brent"/> sur un sous-intervalle autour du meilleur point.
	/// </remarks>
	[[nodiscard]] double Estimate(const astig_model& model, const Vec2D& sigma, const std::array<double, 2>& zLimit = {-800.0, 800.0}, size_t zStep = 201) const;

protected:
	double pixelSize_; ///< Taille de pixel en nanomètre. Attention entre les nanomètres et les micromètres lors des appels
	double sampling_;  ///< Techniquement inutile, on n'enregistre pas les coordonnées apres upsampling donc il est toujours à 1.

	//==================== Hyperparamètres de la calibration ====================
	static constexpr size_t N_PARAMS = 5,     ///< Nombre de paramètres par axe.
	                        N_ITER   = 200;   ///< Itérations max LM (par axe).
	static constexpr double LAMBDA_0 = 1e-3,  ///< Lambda initial LM.
	                        TOL_STEP = 1e-9,  ///< Tolérance numérique générique (évite division par 0).
	                        TOL_COST = 1e-12; ///< Critère d'arrêt relatif sur le coût.

	static constexpr size_t BRENT_ITER = 80;   ///< Itérations max Brent.
	static constexpr double BRENT_TOL  = 1e-3, ///< Tolérance (nm) sur l'argmin.
	                        Z_TOL      = 1e-9; ///< Tolérance (nm) pour tests d'égalité flottante sur Z.

	//==================== Modèle / coût ====================
	/// <summary>Calcule sigma (px) pour un axe à une position z (nm).\n
	///
	/// \f[
	///   \text{Avec} \quad u=\dfrac{z-z0}{w}, \quad
	///   \sigma_{nm}=A\sqrt{1+u^2+c_3u^3+c_4u^4},\quad
	///   \sigma_{px}=\dfrac{sampling}{pixelSize}\sigma_{nm}
	/// \f]
	///
	/// </summary>
	/// <param name="coeffs">Coefficients d'astigmatisme pour l'axe considéré : <c>[z0, w, c3, c4, A]</c>.</param>
	/// <param name="z">Position axiale en nanomètre.</param>
	/// <returns>Sigma prédit (pixels).</returns>
	/// <remarks>
	/// Sécurités numériques :
	/// - <c>w</c> doit rester > 0 (sinon singularité).
	/// - Le polynôme peut devenir négatif : on clampe à 0 avant la racine.
	/// </remarks>
	[[nodiscard]] double getSigma(const std::array<double, 5>& coeffs, double z) const;

	/// <summary>Coût pour la minimisation: erreur quadratique \f$\sigma_x/\sigma_y\f$ (\f${px}^2\f$).</summary>
	/// <param name="model">Coefficients d'astigmatisme pour les deux axes.</param>
	/// <param name="z">Position axiale testée en nanomètre.</param>
	/// <param name="sigma">Écart-type mesuré en X/Y (pixels).</param>
	/// <returns>Erreur calculée : \f$(\Delta\sigma_x)^2 + (\Delta\sigma_y)^2\f$.</returns>
	[[nodiscard]] double getError(const astig_model& model, double z, const Vec2D& sigma) const;

	//==================== Calibration LM ====================
	/// <summary> Ajuste un axe (X ou Y) via Levenberg–Marquardt en minimisant la SSE : \f[ \sum_i (\sigma_{pred}(z_i; p) - \sigma_i)^2 \f]</summary>
	/// <param name="sigma">Mesures \f$\sigma\f$ (pixels), taille N.</param>
	/// <param name="z">Positions z (nm), taille N.</param>
	/// <returns>Paramètres ajustés <c>[z0, w, c3, c4, A]</c>.</returns>
	/// <remarks>
	/// - Initialise p par heuristiques (\f$\sigma_{min} \rightarrow z_0\f$, \f$\Delta{z} \rightarrow w\f$, etc.).
	/// - Contraintes : <c>w</c> et <c>A</c> clampés à un epsilon strictement positif.
	/// </remarks>
	[[nodiscard]] std::array<double, 5> fitAxis(const std::vector<double>& sigma, const std::vector<double>& z) const;

	/// <summary>Évalue résidus et jacobien (analytique) pour un axe, et retourne la SSE.</summary>
	/// <param name="sigma">Mesures \f$\sigma_x\f$ (px).</param>
	/// <param name="z">Positions z (nm).</param>
	/// <param name="p">Paramètres courants.</param>
	/// <param name="r">Sortie : résidus \f$ r_i = \sigma_{pred} - \sigma \f$.</param>
	/// <param name="j">Sortie : jacobien \f$ j_{ik} = \dfrac{\partial\sigma_{pred}}{\partial p_k}\f$.</param>
	/// <returns>SSE (\f${px}^2\f$).</returns>
	[[nodiscard]] double lm(const std::vector<double>& sigma, const std::vector<double>& z,
	                        const std::array<double, N_PARAMS>& p, std::vector<double>& r, std::vector<std::array<double, N_PARAMS>>& j) const;

	/// <summary>Calcule la SSE (somme des carrés des résidus) pour un axe et un jeu de paramètres.</summary>
	/// <param name="sigma">Mesures \f$\sigma_x\f$ (px).</param>
	/// <param name="z">Positions z (nm).</param>
	/// <param name="p">Paramètres courants.</param>
	/// <returns>SSE (\f${px}^2\f$).</returns>
	[[nodiscard]] double computeSSE(const std::vector<double>& sigma, const std::vector<double>& z, const std::array<double, N_PARAMS>& p) const;

	/// <summary>Résout le système linéaire \f$A x = b\f$ (\f$5 \times 5\f$) par élimination de Gauss avec pivot partiel.</summary>
	/// <param name="a">Matrice A (modifiée en place).</param>
	/// <param name="b">Second membre (modifié en place).</param>
	/// <param name="x">Solution.</param>
	/// <returns><c>true</c> si succès, <c>false</c> si matrice singulière (pivot ~ 0).</returns>
	static bool solver(std::array<std::array<double, N_PARAMS>, N_PARAMS>& a, std::array<double, N_PARAMS>& b, std::array<double, N_PARAMS>& x);

	//==================== Estimation Z ====================
	/// <summary>Recherche du minimum de <c>getError(model, z, sigma)</c> sur l'intervalle [a, b] (méthode de Brent).</summary>
	/// <remarks>
	/// La méthode de Brent combine :
	/// - une étape "rapide" : interpolation parabolique (quand la fonction est suffisamment régulière)
	/// - une étape "robuste" : section dorée (garantie de convergence même si la parabole est mauvaise)
	///
	/// Ici, la fonction minimisée est l'erreur entre les sigmas prédits par le modèle et les sigmas mesurés.
	///
	/// Hypothèses / limites importantes :
	/// - La fonction n'a pas besoin d'être strictement convexe sur [a, b].
	/// - Il est toutefois préférable d'appeler Brent sur un intervalle qui contient un unique minimum "local" (d'où l'intérêt d'un pré-scan avant Brent dans Estimate()).
	/// </remarks>
	/// <param name="a">Borne basse (Z en nm).</param>
	/// <param name="b">Borne haute (Z en nm).</param>
	/// <param name="model">Modèle d'astigmatisme (2×5 coefficients).</param>
	/// <param name="sigma">Sigmas mesurés (X/Y en px).</param>
	/// <returns>La valeur de z (nm) qui minimise l'erreur sur [a, b].</returns>
	[[nodiscard]] double brent(double a, double b, const astig_model& model, const Vec2D& sigma) const;
};
} // namespace palm_tracer::fit
