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
using astig_model = std::array<std::array<double, 5>, 2>;

/// <summary>Ajustement non linéaire (Levenberg–Marquardt) d'une courbe d'astigmatisme \f$\sigma(z)\f$, et estimation de Z à partir d'une paire \f$(\sigma_x,
/// \sigma_y) \f$.</summary>
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
	/// <param name="center">Permet de centrer le modèle pour que \f$\sigma_x(0) \approx \sigma_y(0)\f$.</param>
	/// <returns>
	/// Le modèle 2×5 :
	/// - <c>res[0] = [z0x, wx, c3x, c4x, Ax]</c>
	/// - <c>res[1] = [z0y, wy, c3y, c4y, Ay]</c>
	/// </returns>
	/// <remarks>
	/// - Les axes X et Y sont ajustés indépendamment (<see cref="fitAxis"/>).
	/// - En cas d'entrée invalide (tailles incohérentes, N trop petit), l'axe concerné renvoie un vecteur rempli de -1.
	/// - Les candidats sont contrôlés sur l'intervalle des z fournis, avant recentrage. Un recentrage de Zc déplace
	/// aussi cet intervalle vers [min(z)-Zc, max(z)-Zc] ; il ne garantit pas une extrapolation jusqu'à un Zmax arbitraire.
	/// </remarks>
	[[nodiscard]] astig_model Calibrate(const std::vector<Vec3D>& input, bool center = true) const;

	/// <summary>Estime z (nm) en minimisant l'erreur quadratique entre \f$(\sigma_x, \sigma_y)\f$ mesurés et prédits par le modèle.</summary>
	/// <param name="model">Modèle d'astigmatisme (résultat de <see cref="Calibrate"/>).</param>
	/// <param name="sigma">Mesure \f$(\sigma_x, \sigma_y)\f$ en pixels.</param>
	/// <param name="zLimit">Intervalle de recherche initial : <c>[zMin, zMax]</c> en nm.</param>
	/// <param name="zStep">Nombre d'échantillons pour le pré-scan (grille grossière) avant raffinement.
	/// Doit être ≥ 3 (sinon division par <c>zStep-1</c>) et impaire si l'on veut le point 0 dans les choix initiaux.</param>
	/// <returns><c>x</c> : Z (nm) minimisant la SSE ; <c>y</c> : score de compatibilité des sigmas entre 0 et 1, plus grand = meilleur.</returns>
	/// <remarks>
	/// Stratégie :
	/// 1) pré-scan sur <paramref name="zLimit"/> pour trouver un bon point initial,
	/// 2) raffinement par <see cref="brent"/> sur un sous-intervalle autour du meilleur point.
	/// Le score est calculé ensuite par <see cref="getQuality"/> ; il ne participe pas à la minimisation.
	/// Il remplace l'ancienne SSE retournée, sans changer la position Z. Voir \ref axial_quality_score pour la formule et ses limites.
	/// </remarks>
	[[nodiscard]] Vec2D Estimate(const astig_model& model, const Vec2D& sigma,
	                             const std::array<double, 2>& zLimit = {-800.0, 800.0}, size_t zStep = 201) const;

	/// <summary>Estime z (nm) en minimisant l'erreur quadratique entre \f$(\sigma_x, \sigma_y)\f$ mesurés et prédits par le modèle.</summary>
	/// <param name="model">Modèle d'astigmatisme (résultat de <see cref="Calibrate"/>).</param>
	/// <param name="sigmas">Vecteur de mesures \f$(\sigma_x, \sigma_y)\f$ en pixels.</param>
	/// <param name="zLimit">Intervalle de recherche initial : <c>[zMin, zMax]</c> en nm.</param>
	/// <param name="zStep">Nombre d'échantillons pour le pré-scan (grille grossière) avant raffinement.
	/// Doit être ≥ 3 (sinon division par <c>zStep-1</c>) et impaire si l'on veut le point 0 dans les choix initiaux.</param>
	/// <returns>Pour chaque mesure, <c>x</c> = Z (nm) et <c>y</c> = score de compatibilité entre 0 et 1, comme la surcharge scalaire.</returns>
	/// <remarks>
	/// Stratégie :
	/// 1) pré-scan sur <paramref name="zLimit"/> pour trouver un bon point initial,
	/// 2) raffinement par <see cref="brent"/> sur un sous-intervalle autour du meilleur point.
	/// </remarks>
	[[nodiscard]] std::vector<Vec2D> Estimate(const astig_model& model, const std::vector<Vec2D>& sigmas,
	                                          const std::array<double, 2>& zLimit = {-800.0, 800.0}, size_t zStep = 201) const;

protected:
	double pixelSize_; ///< Taille de pixel en nanomètre. Attention entre les nanomètres et les micromètres lors des appels
	double sampling_;  ///< Techniquement inutile, on n'enregistre pas les coordonnées après upsampling donc il est toujours à 1.

	//==================== Hyperparamètres de la calibration ====================
	static constexpr size_t N_PARAMS = 5,     ///< Nombre de paramètres par axe.
	                        N_ITER   = 200;   ///< Itérations max LM (par axe).
	static constexpr double LAMBDA_0 = 1e-3,  ///< Lambda initial LM.
	                        TOL_STEP = 1e-9,  ///< Tolérance numérique générique (évite division par 0).
	                        TOL_COST = 1e-12; ///< Critère d'arrêt relatif sur le coût.

	static constexpr size_t BRENT_ITER    = 80;   ///< Itérations max Brent.
	static constexpr double BRENT_TOL     = 1e-3, ///< Tolérance (nm) sur l'argmin.
	                        Z_TOL         = 1e-9; ///< Tolérance (nm) pour tests d'égalité flottante sur Z.
	static constexpr double QUALITY_SCALE = 0.2;  ///< Écart relatif RMS de 20 % correspondant au score 0.5 ; convention, pas une incertitude.

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

	/// <summary>Convertit les résidus relatifs des sigmas en un score local de compatibilité, sans modifier le coût d'optimisation.
	/// \f[r^2=\frac12\sum_{a\in\{x,y\}}\left(\frac{\sigma_a-\hat\sigma_a(z)}{\hat\sigma_a(z)}\right)^2,
	/// \qquad Q_Z=\frac{1}{1+r^2/\tau^2},\quad\tau=0.20.\f]</summary>
	/// <param name="model">Modèle astigmatique dont sont issues les largeurs prédites.</param>
	/// <param name="z">Position Z retenue par la minimisation de la SSE, en nm.</param>
	/// <param name="sigma">Largeurs observées, dans la même unité que les prédictions (pixels).</param>
	/// <returns>Score entre 0 et 1 ; 0 si une largeur observée ou prédite est non finie ou non positive.</returns>
	/// <remarks>Convention propre au projet, distincte d'un R2, d'une probabilité et d'une précision axiale.
	/// Ne dépend ni des autres localisations ni de la ROI ; <c>zLimit</c> intervient seulement dans le choix de Z.
	/// Un score élevé ne garantit ni une solution unique ni une bonne sensibilité axiale. Voir \ref axial_quality_score.
	/// </remarks>
	[[nodiscard]] double getQuality(const astig_model& model, double z, const Vec2D& sigma) const;

	//==================== Calibration LM ====================
	/// <summary> Ajuste un axe (X ou Y) via Levenberg–Marquardt en minimisant la SSE : \f[ \sum_i (\sigma_{pred}(z_i; p) - \sigma_i)^2 \f]</summary>
	/// <param name="sigma">Mesures \f$\sigma\f$ (pixels), taille N.</param>
	/// <param name="z">Positions z (nm), taille N.</param>
	/// <returns>Paramètres ajustés <c>[z0, w, c3, c4, A]</c>.</returns>
	/// <remarks>
	/// - Initialise p par heuristiques (\f$\sigma_{min} \rightarrow z_0\f$, \f$\Delta{z} \rightarrow w\f$, etc.).
	/// - Contraintes : <c>w</c> et <c>A</c> clampés à un epsilon strictement positif.
	/// - Un candidat n'est accepté que si son polynôme reste strictement supérieur à 1e-18 sur tout
	/// l'intervalle de calibration [min(z), max(z)], afin de rester hors de la troncature du Jacobien.
	/// Cette contrainte ne garantit pas la validité d'une extrapolation au-delà de cet intervalle.
	/// </remarks>
	[[nodiscard]] std::array<double, 5> fitAxis(const std::vector<double>& sigma, const std::vector<double>& z) const;

	/// <summary>Vérifie la positivité du polynôme d'un axe sur un intervalle fermé de calibration.</summary>
	/// <param name="p">Paramètres [z0, w, c3, c4, A], finis, avec w et A strictement positifs.</param>
	/// <param name="zMin">Borne inférieure de calibration, en nm.</param>
	/// <param name="zMax">Borne supérieure de calibration, en nm.</param>
	/// <returns>True si le polynôme est fini et supérieur à 1e-18 aux bornes et à tous ses points stationnaires intérieurs.</returns>
	/// <remarks>Avec \f$P(u)=1+u^2+c_3u^3+c_4u^4\f$, les points stationnaires vérifient
	/// \f$P'(u)=u(2+3c_3u+4c_4u^2)=0\f$. Le contrôle est donc continu, sans grille ni contrainte de symétrie.</remarks>
	[[nodiscard]] static bool isValidAxis(const std::array<double, N_PARAMS>& p, double zMin, double zMax);

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
	/// <param name="a">Borne basse Z (nm).</param>
	/// <param name="b">Borne haute Z (nm).</param>
	/// <param name="model">Modèle d'astigmatisme (2×5 coefficients).</param>
	/// <param name="sigma">Sigmas mesurés (X/Y en px).</param>
	/// <returns>La valeur de Z (nm) qui minimise l'erreur sur [a, b] ainsi que l'erreur associée.</returns>
	[[nodiscard]] Vec2D brent(double a, double b, const astig_model& model, const Vec2D& sigma) const;

	/// <summary>Recherche la position axiale où les courbes SigmaX et SigmaY se croisent le mieux.</summary>
	/// <param name="model">Modèle astigmatique complet [X, Y].</param>
	/// <param name="z">Positions axiales expérimentales en nm.</param>
	/// <returns>
	/// Position axiale du centre astigmatique, en nm.
	/// Si aucun croisement strict n'est trouvé, renvoie la position minimisant |SigmaX - SigmaY|.
	/// </returns>
	[[nodiscard]] double findAstigmaticCenter(const astig_model& model, const std::vector<double>& z) const;
};
} // namespace palm_tracer::fit
