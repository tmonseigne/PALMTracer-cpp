///-------------------------------------------------------------------------------------------------
///
/// \file Fitting.hpp
/// \brief Classe mère pour l'ajustement (fitting) de PSF sur ROI 2D/3D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 27/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <cstdint>
#include <vector>

#include "Constants.hpp"
#include "Structs.hpp"
#include "localization/Segmentation2D.hpp"

namespace palm_tracer::fit {
/// <summary>Mode d'ajustement.</summary>
enum class fit_mode : uint8_t {
	none,              ///< Aucun ajustement.
	gaussian_2d_fixed, ///< Ajustement Gaussien 2D avec \f$\sigma\f$ fixé (uniquement X et Y).
	gaussian_2d_sigma, ///< Ajustement Gaussien 2D avec estimation de \f$\sigma\f$ isotrope.
	gaussian_3d_fixed, ///< Ajustement Gaussien 3D avec \f$\theta\f$ fixé et estimation de \f$\sigma_x\f$ et \f$\sigma_y\f$.
	gaussian_3d_theta, ///< Ajustement Gaussien 3D avec estimation de \f$\sigma_x\f$, \f$\sigma_y\f$ et \f$\theta\f$</c>.
	spline_3d_emccd,   ///< Ajustement par PSF spline (capteur EMCCD).
	spline_3d_scmos    ///< Ajustement par PSF spline (capteur sCMOS).
};

/// <summary>Classe d'ajustement gaussien/spline générique.</summary>
class Fitting {
public:
	/// <summary>Constructeur de la classe <see cref="Fitting"/>.\n
	/// Initialise le nombre de paramètres du modèle en fonction du mode d'ajustement choisi.</summary>
	/// <param name="mode">Mode de fitting gaussien/spline.</param>
	/// <param name="image">Image 2D d'entrée (hauteur × largeur), valeurs entières non-signées 16 bits.</param>
	/// <param name="dim">Dimensions de l'image (X, Y).</param>
	/// <param name="maxPoints">Nombre maximum de points par plan (prévention d'un débordement dans la sortie <c>points</c> de <see cref="Fit"/>).</param>
	/// <param name="fitParams">Paramètres pour l'ajustement (dépendent du <paramref name="mode"/> ; voir <see cref="parseParams"/>).</param>
	explicit Fitting(const fit_mode mode, const uint16_t* image, const Vec2I& dim, const size_t maxPoints, const double* fitParams)
		: mode_(mode), image_(image), dim_(dim), maxPoints_(maxPoints), fitParams_(fitParams) { initialize(); }

	/// <summary>Destructeur par défaut de la classe <see cref="Fitting"/>.</summary>
	virtual ~Fitting() = default;

	/// <summary>Exécute l'ajustement sur tous les points préalablement localisés par la segmentation.</summary>
	/// <param name="segmentation">Segmentation de l'image initiale (liste d'objets à ajuster).</param>
	/// <param name="points">Sortie : tableau des localisations ajustées (taille ajustée par la méthode).</param>
	/// <remarks>
	/// Parcourt les objets retournés par <see cref="Segmentation2D::Infos"/> et applique <see cref="Fit"/> à chacun.
	/// Bouclage parallèle (<c>OpenMP</c>) jusqu'à <see cref="maxPoints_"/>.
	/// Émet un avertissement sur la sortie d'erreur une seule fois par appel si des points dépassent la capacité.
	/// </remarks>
	void FitAll(const localization::Segmentation2D& segmentation, std::vector<std::array<double, N_LOCALIZATIONS>>& points) const;

	/// <summary>Exécute l'ajustement sur une seule région, sans créer d'équipe OpenMP.</summary>
	/// <param name="infos">Informations issues de la segmentation (position initiale, intensité, etc.).</param>
	/// <param name="point">Sortie : localisation ajustée (format <see cref="N_LOCALIZATIONS"/>).</param>
	/// <remarks>Implémentation par défaut : recopie les valeurs brutes via <see cref="fill"/>.
	/// L'appelant peut répartir les régions d'un lot entre ses propres threads.</remarks>
	virtual void Fit(const localization::SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) const
	{
		fill(infos, point, infos.x_, infos.y_, infos.z_, 0, 0);
	}

protected:
	fit_mode mode_;           ///< Mode de fitting utilisé.
	const uint16_t* image_;   ///< Image 2D d'entrée.
	Vec2I dim_;               ///< Dimensions de l'image (X, Y).
	size_t n_{0},             ///< Nombre d'échantillons (= \f$roi \times roi\f$).
	       nParams_{4},       ///< Nombre total de paramètres du modèle.
	       nParams2_{},       ///< Nombre de paramètres au carré (dimension de matrices de type \f$J^\top J\f$).
	       maxPoints_,        ///< Nombre maximal de points traités par plan (sécurité).
	       roi_{},            ///< Taille de la région d'intérêt autour du point (toujours impaire, centrée).
	       roiHalf_{};        ///< Demi-taille de la région d'intérêt.
	const double* fitParams_; ///< Paramètres utilisateur pour le fitting.

	static constexpr size_t N_ITER = 75; ///< Nombre maximal d'itérations autorisées par l'optimiseur.

	/// <summary>Initialise les membres et valide les paramètres.
	/// Rend <see cref="roi_"/> impair (si pair, on décrémente), met à jour <see cref="roiHalf_"/>, <see cref="n_"/> et <see cref="nParams2_"/>.</summary>
	void initialize();

	/// <summary>Analyse les paramètres utilisateur (<see cref="fitParams_"/>) et configure l'objet (<see cref="roi_"/>, etc.).</summary>
	/// <remarks>Implémentation par défaut : <c>roi_ = static_cast&lt;size_t&gt;(fitParams_[0])</c>. Les classes dérivées peuvent étendre cette logique.</remarks>
	virtual void parseParams() { roi_ = static_cast<size_t>(fitParams_[0]); }

	/// <summary>Remplit le vecteur de sortie <paramref name="point"/> après un ajustement.</summary>
	/// <param name="infos">Informations initiales sur le point localisé.</param>
	/// <param name="point">Sortie : localisation au format interne.</param>
	/// <param name="x">Coordonnée X estimée (en pixels).</param>
	/// <param name="y">Coordonnée Y estimée (en pixels).</param>
	/// <param name="z">Coordonnée Z estimée (unités dépendant du modèle spline ; 0 par défaut).</param>
	/// <param name="offset">Fond/offset estimé.</param>
	/// <param name="intensity">Amplitude/Intensité estimée.</param>
	/// <param name="sigmaX">Écart-type estimé sur X (pixels).</param>
	/// <param name="sigmaY">Écart-type estimé sur Y (pixels).</param>
	/// <param name="theta">Angle d'ellipse (degrés, convention image).</param>
	/// <param name="rssXY">Erreur quadratique moyenne (MSE) sur X et Y.</param>
	/// <param name="rssZ">Erreur quadratique moyenne (MSE) sur Z (si applicable).</param>
	void fill(const localization::SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point, double x, double y, double z = 0,
	          double offset = -1, double intensity = -1, double sigmaX = -1, double sigmaY = -1, double theta = -1,
	          double rssXY = -1, double rssZ = -1) const;

	/// <summary>Prépare les données de la ROI centrée en (<paramref name="x"/>, <paramref name="y"/>) :
	/// coordonnées relatives, intensités et bornes [min, max].</summary>
	/// <param name="x">Coordonnée X centrale de la ROI (index colonne).</param>
	/// <param name="y">Coordonnée Y centrale de la ROI (index ligne).</param>
	/// <param name="grid">Sortie : grille des pixels (<c>x, y, valeur</c>, <see cref="Pixel"/>).</param>
	/// <param name="min">Sortie : intensité minimale détectée dans la ROI.</param>
	/// <param name="max">Sortie : intensité maximale détectée dans la ROI.</param>
	void preCompute(size_t x, size_t y, std::vector<Pixel>& grid, double& min, double& max) const;

	//==================== Run Computes ====================
	/// <summary>Initialise les paramètres de départ (<c>p[1]</c>) ainsi que leurs bornes min/max (<c>p[0]</c>/<c>p[2]</c>) pour l'optimisation.
	/// Le contenu dépend du <see cref="fit_mode"/> courant.</summary>
	/// <param name="p">Tableau de 3 vecteurs : <c>p[0]=min</c>, <c>p[1]=init</c>, <c>p[2]=max</c>.</param>
	/// <param name="x">Estimation initiale de X.</param>
	/// <param name="y">Estimation initiale de Y.</param>
	/// <param name="z">Estimation initiale de Z (modèles spline).</param>
	/// <param name="minI">Intensité minimale mesurée dans la ROI (sert à initialiser l'offset).</param>
	/// <param name="maxI">Intensité maximale mesurée dans la ROI (sert à initialiser l'amplitude).</param>
	void initP(std::array<std::vector<double>, 3>& p, double x, double y, double z, double minI, double maxI) const;
};
} // namespace palm_tracer::fit
