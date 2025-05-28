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
	explicit GaussianFit(const fit_mode mode, const unsigned short* image, const size_t height, const size_t width, const size_t maxPoints,
	                     const size_t roi, const double sigmaX, const double sigmaY, const double theta)
		: Fitting(mode, image, height, width, maxPoints, roi, sigmaX, sigmaY, theta) { }

protected:
	/// <summary>Exécute un ajustement de modèle Gaussien sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.
	/// L'algorithme utilise une méthode de Levenberg-Marquardt pour minimiser l'erreur entre le modèle et les données.</summary>
	void fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point) override;


	//==================== Run Computes ====================
	/// <summary>Prépare les données nécessaires pour le fitting en extrayant les valeurs de la région d'intérêt autour de (x, y).
	/// Remplit les vecteurs de données, poids et coordonnées relatives. Calcule aussi les intensités minimale et maximale.</summary>
	/// <param name="min">Référence vers la valeur minimale trouvée dans la ROI (définie par la méthode).</param>
	/// <param name="max">Référence vers la valeur maximale trouvée dans la ROI (définie par la méthode).</param>
	/// <param name="image">Pointeur vers les données de l'image source (unsigned short).</param>
	/// <param name="width">Largeur de l'image (nombre de colonnes).</param>
	/// <param name="x">Coordonnée X centrale de la ROI.</param>
	/// <param name="y">Coordonnée Y centrale de la ROI.</param>
	/// <param name="dat">Vecteur de sortie contenant les données d'intensité extraites de la ROI.</param>
	/// <param name="weight">Vecteur de poids initialisés pour chaque pixel de la ROI (initialement tous à 1).</param>
	/// <param name="t">Vecteur contenant les coordonnées X et Y associées à chaque point de la ROI.</param>
	/// <param name="roi">Taille (impair) du carré de la région d'intérêt à analyser.</param>
	static void preComputeFit(double& min, double& max, const unsigned short* image, size_t width, size_t x, size_t y,
	                          std::vector<double>& dat, std::vector<double>& weight, std::vector<size_t>& t, size_t roi);

	/// <summary>Initialise les paramètres de départ pour le fitting, ainsi que leurs bornes minimale et maximale.
	/// Les valeurs dépendent du mode de fitting utilisé.</summary>
	/// <param name="p"></param>
	/// <param name="x">Coordonnée X estimée du centre initial.</param>
	/// <param name="y">Coordonnée Y estimée du centre initial.</param>
	/// <param name="sigmaX">Valeur initiale de la largeur gaussienne en X.</param>
	/// <param name="sigmaY">Valeur initiale de la largeur gaussienne en Y.</param>
	/// <param name="theta">Valeur initiale de la rotation (en radians, utilisée uniquement dans le mode 3D avec angle).</param>
	/// <param name="minI">Intensité minimale dans la ROI (servira à initialiser l'offset).</param>
	/// <param name="maxI">Intensité maximale dans la ROI (servira à initialiser l'amplitude).</param>
	void initP(std::array<std::vector<double>, 3>& p, double x, double y, double sigmaX, double sigmaY, double theta, double minI, double maxI) const;

	/// <summary>Calcule le coefficient de corrélation au carré (R²) entre deux vecteurs normalisés.
	/// Cela donne une mesure de la qualité d'ajustement entre les données observées et le modèle.</summary>
	/// <param name="dat">Vecteur des données observées (modifié en place).</param>
	/// <param name="hat">Vecteur des données ajustées (modifié en place).</param>
	/// <returns>Valeur de R², dans [0, 1].</returns>
	double computeR2(std::vector<double>& dat, std::vector<double>& hat) const;

	//==================== PSF ====================
	/// <summary>Évalue la fonction PSF (Point Spread Function) en (x, y), selon le mode d'ajustement et les paramètres donnés.</summary>
	/// <param name="p">Vecteur des paramètres de la PSF (intensité, position, sigma, offset...).</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	/// <returns>Valeur de la PSF en (x, y).</returns>
	double psf(const std::vector<double>& p, double x, double y) const;

	/// <summary>PSF gaussienne isotrope (2D) centrée en (0,0), à sigma fixe.</summary>
	/// <param name="x">Coordonnée x relative au centre.</param>
	/// <param name="y">Coordonnée y relative au centre.</param>
	/// <param name="i0">Amplitude de la gaussienne.</param>
	/// <param name="sigma">Écart-type commun (isotrope).</param>
	/// <param name="offset">Valeur de fond (background constant).</param>
	/// <returns>Valeur de la PSF en (x, y).</returns>
	static double psf(double x, double y, double i0, double sigma, double offset);

	/// <summary>PSF gaussienne anisotrope 2D avec rotation (theta).</summary>
	/// <param name="x">Coordonnée x relative au centre.</param>
	/// <param name="y">Coordonnée y relative au centre.</param>
	/// <param name="i0">Amplitude de la gaussienne.</param>
	/// <param name="sigmaX">Écart-type selon l'axe X.</param>
	/// <param name="sigmaY">Écart-type selon l'axe Y.</param>
	/// <param name="theta">Angle de rotation (en radians).</param>
	/// <param name="offset">Valeur de fond.</param>
	/// <returns>Valeur de la PSF en (x, y).</returns>

	static double psf(double x, double y, double i0, double sigmaX, double sigmaY, double theta, double offset);

	//==================== Matrix ====================
	/// <summary>Calcule le gradient de la PSF par rapport aux paramètres, pour un point (x, y), selon le mode d'ajustement.</summary>
	/// <param name="out">Pointeur vers le tableau de sortie du gradient (taille dépendant du mode).</param>
	/// <param name="p">Paramètres de la PSF.</param>
	/// <param name="x">Coordonnée x d'évaluation.</param>
	/// <param name="y">Coordonnée y d'évaluation.</param>
	void deriv(double* out, const std::vector<double>& p, double x, double y) const;

	/// <summary>Dérivées partielles de la PSF 2D isotrope (centrée) par rapport aux paramètres.</summary>
	/// <param name="out">Tableau de sortie (au moins 5 éléments).</param>
	/// <param name="x">Coordonnée x relative au centre.</param>
	/// <param name="y">Coordonnée y relative au centre.</param>
	/// <param name="i0">Amplitude.</param>
	/// <param name="sigma">Écart-type commun.</param>
	/// <param name="has_sigma">True si sigma est un paramètre libre (ajouter dérivée), false sinon.</param>
	static void deriv(double* out, double x, double y, double i0, double sigma, bool has_sigma);

	/// <summary>Dérivées partielles de la PSF 2D anisotrope avec rotation par rapport aux paramètres.</summary>
	/// <param name="out">Tableau de sortie (au moins 7 éléments).</param>
	/// <param name="x">Coordonnée x relative au centre.</param>
	/// <param name="y">Coordonnée y relative au centre.</param>
	/// <param name="i0">Amplitude.</param>
	/// <param name="sigmaX">Écart-type selon X.</param>
	/// <param name="sigmaY">Écart-type selon Y.</param>
	/// <param name="theta">Angle de rotation.</param>
	/// <param name="has_theta">True si theta est un paramètre libre, false sinon.</param>
	static void deriv(double* out, double x, double y, double i0, double sigmaX, double sigmaY, double theta, bool has_theta);

	/// <summary>Calcule les matrices alpha et beta du système linéaire de l'ajustement de Levenberg-Marquardt, ainsi que la valeur de chi² pondéré entre modèle et données.</summary>
	/// <param name="alpha">Matrice des coefficients (nParams x nParams), remplie par la fonction.</param>
	/// <param name="beta">Vecteur des termes constants du système linéaire.</param>
	/// <param name="dydp">Dérivées partielles des valeurs ajustées par rapport aux paramètres.</param>
	/// <param name="weight">Poids associés à chaque point.</param>
	/// <param name="deltaY">Résidus entre données observées et ajustées.</param>
	/// <returns>Valeur du chi² pondéré.</returns>
	double chi2(std::vector<double>& alpha, std::vector<double>& beta, const std::vector<double>& dydp,
	            const std::vector<double>& weight, const std::vector<double>& deltaY) const;

	/// <summary>Calcule la matrice alpha, le vecteur beta, et retourne le chi² (erreur) pour l'itération courante.
	/// Utilisé dans l'algorithme de Levenberg-Marquardt.</summary>
	/// <param name="t">Vecteur de coordonnées (x, y) pour chaque point</param>
	/// <param name="p">Paramètres du modèle</param>
	/// <param name="dat">Données observées</param>
	/// <param name="weight">Poids pour chaque point</param>
	/// <param name="alpha">Matrice approximative d'Hessienne</param>
	/// <param name="beta">Gradient du χ²</param>
	/// <param name="hat">Vecteur des données ajustées (modifié en place).</param>
	/// <returns>Valeur du χ²</returns>
	double lm_matx(const std::vector<size_t>& t, const std::vector<double>& p, const std::vector<double>& dat, const std::vector<double>& weight,
	               std::vector<double>& alpha, std::vector<double>& beta, std::vector<double>& hat) const;

	/// <summary>Inversion de matrice (méthode spécifique à l'implémentation, incomplète ici).</summary>
	/// <param name="in">Matrice à inverser (format vectorisé ligne par ligne).</param>
	/// <param name="out">Matrice résultat de l'inversion (même format).</param>
	void inverse(std::vector<double>& in, std::vector<double>& out) const;
};
}
