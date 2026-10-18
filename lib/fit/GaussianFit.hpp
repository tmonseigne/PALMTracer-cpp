///-------------------------------------------------------------------------------------------------
///
/// \file GaussianFit.hpp
/// \brief Classe d'ajustement gaussien (2D et 3D) par Levenberg–Marquardt sur ROI.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 26/03/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <vector>

#include "Fitting.hpp"

namespace palm_tracer::fit {
/// <summary>Classe d'ajustement Gaussien dérivée de <see cref="Fitting"/>.
/// Elle réalise un ajustement non linéaire (LM) d'un modèle de PSF Gaussienne sur une ROI carrée centrée sur le candidat, pour plusieurs <see cref="fit_mode"/>.\n
/// Modes pris en charge :
/// - <see cref="fit_mode::gaussian_2d_fixed"/> : \f$p=[\text{offset}, I, x, y]\f$, \f$\sigma_x=\sigma_y\f$ fixé.
/// - <see cref="fit_mode::gaussian_2d_sigma"/> : \f$p=[\text{offset}, I, x, y, \sigma]\f$ (isotrope).
/// - <see cref="fit_mode::gaussian_3d_fixed"/> : \f$p=[\text{offset}, I, x, y, \sigma_x\sqrt{2}, \sigma_y\sqrt{2}]\f$, \f$\theta\f$ fixé.
/// - <see cref="fit_mode::gaussian_3d_theta"/> : \f$p=[\text{offset}, I, x, y, \sigma_x\sqrt{2}, \sigma_y\sqrt{2}, \theta]\f$.
///
/// Les \f$\sqrt{2}\f$ proviennent de la convention interne du modèle 3D (cf. normalisation dans <see cref="fit"/>).
/// </summary>
class GaussianFit final : public Fitting {
public:
	/// <summary>Constructeur de la classe <see cref="GaussianFit"/>.</summary>
	/// <param name="mode">Mode d'ajustement (<see cref="fit_mode"/>).</param>
	/// <param name="image">Pointeur vers l'image 16 bits source.</param>
	/// <param name="dim">Dimensions de l'image (en pixels).</param>
	/// <param name="maxPoints">Nombre maximum de points à traiter.</param>
	/// <param name="fitParams">Paramètres de configuration (voir <see cref="parseParams"/>).</param>
	explicit GaussianFit(const fit_mode mode, const uint16_t* image, const Vec2I& dim, const size_t maxPoints,
	                     const double* fitParams) : Fitting(mode, image, dim, maxPoints, fitParams) { initialize(); }

	/// <summary>Ajuste le modèle Gaussien sur la ROI centrée en (<paramref name="infos"/>.x\_, y\_).\n
	/// L'algorithme LM minimise \f$\chi^2=r^\top W r\f$ jusqu'aux critères d'arrêt.
	/// Remplit <paramref name="point"/> avec offset, intensité, (x, y, z), \f$\sigma_x\f$, \f$\sigma_y\f$, \f$\theta\f$, R².</summary>
	/// <param name="infos">Informations issues de la segmentation (position initiale, intensité, etc.).</param>
	/// <param name="point">Sortie : localisation ajustée (format <see cref="N_LOCALIZATIONS"/>).</param>
	void Fit(const localization::SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const override;

protected:
	//==================== Hyperparamètres LM ====================
	static constexpr double LAMBDA_0   = 1e-2, ///< Valeur initiale de \f$\lambda\f$ (LM).
	                        LAMBDA_MIN = 1e-7, ///< Valeur minimale autorisée pour \f$\lambda\f$.
	                        LAMBDA_MAX = 1e7,  ///< Valeur maximale autorisée pour \f$\lambda\f$.
	                        LAMBDA_UP  = 11,   ///< Facteur d'augmentation de \f$\lambda\f$ si l'essai est refusé.
	                        LAMBDA_DN  = 9;    ///< Facteur de diminution de \f$\lambda\f$ si l'essai est accepté.
	static constexpr double EPSILON_1  = 1e-6, ///< Seuil sur \f$\max|\,\beta\,|\f$ (gradient) pour l'arrêt.
	                        EPSILON_2  = 1e-6, ///< Seuil sur la variation relative des paramètres.
	                        EPSILON_3  = 1e-5, ///< Seuil sur l'erreur moyenne par pixel.
	                        EPSILON_4  = 1e-2; ///< Seuil pour accepter un gain (ratio \f$\rho\f$).

	/// <summary>Paramètres initiaux du modèle.</summary>
	double sigmaX_{1.0}, ///< Écart-type initial en X (pixels).
	       sigmaY_{2.0}, ///< Écart-type initial en Y (pixels).
	       theta_{0.0};  ///< Angle d'orientation (radians), selon le mode.

	/// <summary>Espace de travail pour Levenberg–Marquardt.\n
	/// Chaque membre correspond à une variable intermédiaire standard de l'algorithme.</summary>
	struct LM_Workspace {
		std::vector<double> alpha,            ///< Matrice de courbure \f$\alpha=J^\top W J\f$ (taille \f$n_p^2\f$).
		                    alphaR,           ///< Matrice régularisée \f$\alpha_R=\alpha+\lambda\,\mathrm{diag}(\alpha)\f$.
		                    alphaRi,          ///< Inverse de \f$\alpha_R\f$ (Gauss–Jordan).
		                    beta,             ///< Vecteur \f$\beta=J^\top W r\f$ (taille \f$n_p\f$).
		                    delta,            ///< Résidus non pondérés \f$\delta=y-\hat{y}\f$.
		                    deltaW,           ///< Résidus pondérés \f$W\delta\f$.
		                    deltaP,           ///< Incrément \f$\delta p\f$.
		                    lBetaDeltaP,      ///< Terme \f$\lambda\,\delta p+\beta\f$ pour le calcul du gain.
		                    grad,             ///< Gradient de la PSF au pixel courant.
		                    pTry,             ///< Paramètres candidats \f$p+\delta p\f$ (clamper aux bornes).
		                    hat,              ///< Prédiction \f$\hat{y}\f$ (valeurs du modèle sur la ROI).
		                    weight;           ///< Poids par pixel.
		std::vector<Pixel> grid;              ///< Grille de la ROI (coordonnées + intensité).
		std::array<std::vector<double>, 3> p; ///< Bornes inférieures, valeurs courantes et bornes supérieures.

		/// <summary>Alloue/redimensionne les buffers.</summary>
		/// <param name="n">Nombre de pixels dans la ROI.</param>
		/// <param name="np">Nombre de paramètres du modèle.</param>
		/// <param name="np2">Carré du nombre de paramètres.</param>
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
	/// <summary>Parse <paramref name="fitParams_"/> provenant du constructeur et initialise <see cref="roi_"/>, <see cref="sigmaX_"/>, <see cref="sigmaY_"/>, <see cref="theta_"/>.</summary>
	void parseParams() override;

	//==================== Levenberg–Marquardt ====================
	/// <summary>Évalue la PSF sur la ROI, accumule \f$\alpha\f$ et \f$\beta\f$ (Gauss–Newton)
	/// et retourne \f$\chi^2=\sum_i w_i\,\delta_i^2\f$ pour l'état courant \f$p\f$.</summary>
	/// <param name="ws">Espace de travail pour Levenberg–Marquardt.</param>
	/// <returns>Valeur de \f$\chi^2\f$.</returns>
	[[nodiscard]] double lm(LM_Workspace& ws) const;

	/// <summary>Calcule le gain prédit \f$g=2\,\delta p^\top \beta-\delta p^\top \alpha\,\delta p\f$.</summary>
	/// <param name="ws">Espace de travail pour Levenberg–Marquardt.</param>
	/// <returns>Gain estimé.</returns>
	[[nodiscard]] double gain(const LM_Workspace& ws) const;

	/// <summary>Calcule \f$\sum_i w_i\,(y_i-\hat{y}_i(p+\delta p))^2\f$ pour décider d'accepter/refuser la mise à jour.</summary>
	/// <param name="ws">Espace de travail pour Levenberg–Marquardt.</param>
	/// <returns>Erreur calculée.</returns>
	[[nodiscard]] double getError(const LM_Workspace& ws) const;

	/// <summary>Calcule \f$R^2\f$ (corrélation au carré) entre données et modèle sur la ROI.</summary>
	/// <param name="ws">Espace de travail pour Levenberg–Marquardt.</param>
	/// <returns>\f$R^2\in[0,1]\f$ (retourne 0 en cas de NaN).</returns>
	[[nodiscard]] double computeR2(LM_Workspace& ws) const;

	//==================== Run Computes ====================
	/// <summary>PSF Gaussienne évaluée en (\f$x,y\f$) selon <see cref="mode_"/>.</summary>
	/// <param name="p">Vecteur des paramètres de la PSF (intensité, position, sigma, offset...).</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	/// <returns>Valeur de la PSF en (x, y).</returns>
	[[nodiscard]] double psf(const std::vector<double>& p, double x, double y) const;

	/// <summary>Gradient de la PSF par rapport aux paramètres (selon <see cref="mode_"/>).</summary>
	/// <param name="grad">Vecteur de sortie du gradient (taille dépendant du mode).</param>
	/// <param name="p">Paramètres de la PSF.</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	void deriv(std::vector<double>& grad, const std::vector<double>& p, double x, double y) const;
};
} // namespace palm_tracer::fit
