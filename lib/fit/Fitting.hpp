///-------------------------------------------------------------------------------------------------
///
/// \file Fitting.hpp
/// \brief Classe mère pour l'ajustement.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <cstdint>
#include <vector>

#include "Constants.hpp"
#include "Utils.hpp"
#include "palm/Segmentation2D.hpp"

namespace palm_tracer {
/// <summary>Mode d'ajustement.</summary>
enum class fit_mode : uint8_t {
	none,              ///< Aucun Ajustement.
	gaussian_2d_fixed, ///< Ajustement Gaussien 2D avec Sigma fixé (uniquement X et Y).
	gaussian_2d_sigma, ///< Ajustement Gaussien 2D avec calcul du Sigma.
	gaussian_3d_fixed, ///< Ajustement Gaussien 3D avec Theta fixé et calcul du Sigma sur X et Y.
	gaussian_3d_theta, ///< Ajustement Gaussien 3D avec calcul du Sigma sur X et Y et du Theta.
	spline_3d_emccd,   ///< Ajustement par Spline pour Capteur EMCCD
	spline_3d_scmos    ///< Ajustement par Spline pour Capteur sCMOS
};

/// <summary>Classe d'ajustement gaussien.</summary>
class Fitting {
public:
	/// <summary>Constructeur de la classe <c>Fitting</c>. Initialise le nombre de paramètres du modèle en fonction du mode d'ajustement choisi.</summary>
	/// <param name="mode">Mode de fitting gaussien (2D, 2D avec sigma, 3D, 3D avec rotation).</param>
	/// <param name="image">Image 2D d'entrée (hauteur x largeur) en niveaux de gris.</param>
	/// <param name="height">Hauteur de l'image.</param>
	/// <param name="width">Largeur de l'image.</param>
	/// <param name="maxPoints">Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).</param>
	/// <param name="fitParams"> Paramètres pour l'ajustement.</param>
	explicit Fitting(const fit_mode mode, const unsigned short* image, const size_t height, const size_t width, const size_t maxPoints, const double* fitParams)
		: mode_(mode), image_(image), n_(0), nParams_(4), height_(height), width_(width), maxPoints_(maxPoints), fitParams_(fitParams) { initialize(); }

	/// <summary>Destructeur de la classe <c>Fitting</c>.</summary>
	virtual ~Fitting() = default;

	/// <summary>Exécute l'ajustement sur tous les points précédemments localisés par la segmentation.</summary>
	/// <param name="segmentation">Segmentation de l'image initiale. </param>
	/// <param name="points">Résultat de l'ajustement sur les points localisés.</param>
	void Fit(const Segmentation2D& segmentation, std::vector<std::array<double, N_LOCALIZATIONS>>& points);

protected:
	fit_mode mode_;               ///< Mode de fitting utilisé (ex. position seulement, ou position + orientation...).
	const unsigned short* image_; ///< Image 2D d'entrée (hauteur x largeur) en niveaux de gris.
	size_t n_,                    ///< Nombre d'échantillons ou points de données à fitter.
	       nParams_,              ///< Nombre total de paramètres du modèle.
	       nParams2_,             ///< Nombre de paramètres au carré (utile pour les matrices jacobiennes ou hessiennes).
	       height_,               ///< Hauteur de l'image.
	       width_,                ///< Largeur de l'image.
	       maxPoints_,            ///< Nombre maximum de points par plan (pour éviter le débordement dans le tableau `points_`).
	       roi_,                  ///< Taille de la région d'intérêt autour du point.
	       roiHalf_;              ///< Demi taille de la région d'intérêt autour du point.
	const double* fitParams_;     ///< Paramètres de l'ajustement

	static constexpr size_t N_ITER = 50; ///< Nombre maximal d'itérations autorisées dans l'algorithme d'optimisation.

	/// <summary>Initialise les différents membres de la classe. </summary>
	void initialize();

	/// <summary>Parse les paramètres reçus pour l'ajustement.</summary>
	virtual void parseParams() { roi_ = 7; }

	/// <summary>Exécute l'ajustement selon le mode sélectionné sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.
	/// L'algorithme utilise une méthode de Levenberg-Marquardt pour minimiser l'erreur entre le modèle et les données.</summary>
	/// <param name="infos">Informations sur le point localisé.</param>
	/// <param name="point">Résultat de l'ajustement sur le point localisé.</param>
	virtual void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) { fill(infos, point, 0, 0, infos.x_, infos.y_, infos.z_); }

	/// <summary></summary>
	/// <param name="infos">Informations initiales sur le point localisé.</param>
	/// <param name="point">Résultat de l'ajustement sur le point localisé.</param>
	/// <param name="offset">Nouvel offset.</param>
	/// <param name="intensity">Nouvelle intensité.</param>
	/// <param name="x">Nouveau X.</param>
	/// <param name="y">Nouveau Y.</param>
	/// <param name="z">Nouveau Z.</param>
	/// <param name="sigmaX">Nouveau Sigma sur X.</param>
	/// <param name="sigmaY">Nouveau Sigma sur Y.</param>
	/// <param name="theta">Nouveau Theta.</param>
	/// <param name="rss">RSS sur X et Y.</param>
	/// <param name="rssZ">RSS sur Z.</param>
	void fill(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point, double offset, double intensity, double x, double y, double z = 0,
	          double sigmaX = 0, double sigmaY = 0, double theta = 0, double rss = -1, double rssZ = -1) const;

	/// <summary>Prépare les données nécessaires pour le fitting en extrayant les valeurs de la région d'intérêt autour de (x, y).
	/// Remplit les vecteurs de données, coordonnées relatives et calcule aussi les intensités minimale et maximale.</summary>
	/// <param name="x">Coordonnée X centrale de la ROI.</param>
	/// <param name="y">Coordonnée Y centrale de la ROI.</param>
	/// <param name="grid">Vecteur contenant les coordonnées X, Y et intensité associée à chaque point de la ROI.</param>
	/// <param name="min">Référence vers la valeur minimale trouvée dans la ROI (définie par la méthode).</param>
	/// <param name="max">Référence vers la valeur maximale trouvée dans la ROI (définie par la méthode).</param>
	void preCompute(size_t x, size_t y, std::vector<Pixel>& grid, double& min, double& max) const;

	//==================== Run Computes ====================
	/// <summary>Initialise les paramètres de départ pour le fitting, ainsi que leurs bornes minimale et maximale.
	/// Les valeurs dépendent du mode de fitting utilisé.</summary>
	/// <param name="p"></param>
	/// <param name="x">Coordonnée X estimée du centre initial.</param>
	/// <param name="y">Coordonnée Y estimée du centre initial.</param>
	/// <param name="z">Coordonnée Z estimée du centre initial.</param>
	/// <param name="minI">Intensité minimale dans la ROI (servira à initialiser l'offset).</param>
	/// <param name="maxI">Intensité maximale dans la ROI (servira à initialiser l'amplitude).</param>
	void initP(std::array<std::vector<double>, 3>& p, double x, double y, double z, double minI, double maxI) const;
};
} // namespace palm_tracer
