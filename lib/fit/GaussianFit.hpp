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
#include <cstdint>
#include <vector>

/// <summary>Mode d'ajustement gaussien.</summary>
enum class fit_mode : uint8_t { gauss_2d, gauss_2d_sigma, gauss_3d, gauss_3d_theta };

/// <summary>Classe d'ajustement gaussien.</summary>
class GaussianFit
{
public:
	/// <summary>Constructeur de la classe <c>GaussianFit</c>. Initialise le nombre de paramètres du modèle en fonction du mode d'ajustement choisi.</summary>
	/// <param name="mode">Mode de fitting gaussien (2D, 2D avec sigma, 3D, 3D avec rotation).</param>
	explicit GaussianFit(fit_mode mode = fit_mode::gauss_2d);

	/// <summary>Destructeur de la classe <c>GaussianFit</c>.</summary>
	~GaussianFit() = default;

	/// <summary>Exécute un ajustement de modèle Gaussien sur une région d'intérêt (ROI) centrée sur une position (x, y) dans l'image fournie.
	/// L'algorithme utilise une méthode de Levenberg-Marquardt pour minimiser l'erreur entre le modèle et les données.</summary>
	/// <param name="image">Pointeur vers les données brutes de l'image (niveau de gris, type unsigned short).</param>
	/// <param name="sizeX">Largeur totale de l'image.</param>
	/// <param name="sizeY">Hauteur totale de l'image.</param>
	/// <param name="x">Coordonnée X du centre du spot à fitter (modifié en sortie avec la valeur ajustée).</param>
	/// <param name="y">Coordonnée Y du centre du spot à fitter (modifié en sortie avec la valeur ajustée).</param>
	/// <param name="intensity">Intensité estimée du spot (valeur crête ajustée, modifiée en sortie).</param>
	/// <param name="sigmaX">Largeur gaussienne en X (modifiable en fonction du mode).</param>
	/// <param name="sigmaY">Largeur gaussienne en Y (modifiable en fonction du mode).</param>
	/// <param name="theta">Angle de rotation du spot (modifié uniquement si mode avec rotation).</param>
	/// <param name="offset">Niveau de fond (offset constant ajouté à la gaussienne).</param>
	/// <param name="rss">Résidu quadratique de l'ajustement (Root Sum of Squares).</param>
	/// <param name="roi">Taille de la région d'intérêt carrée centrée autour de (x, y) utilisée pour le fitting. Doit être impair.</param>
	void Run(const unsigned short* image, size_t sizeX, size_t sizeY, double& x, double& y, double& intensity,
	         double& sigmaX, double& sigmaY, double& theta, double& offset, double& rss, size_t roi);

protected:
	fit_mode mode_;                                    ///< Mode de fitting utilisé (ex. position seulement, ou position + orientation...).
	size_t n_,                                         ///< Nombre d'échantillons ou points de données à fitter.
	       nParams_,                                   ///< Nombre total de paramètres du modèle.
	       nParams2_;                                  ///< Nombre de paramètres au carré (utile pour les matrices jacobiennes ou hessiennes).
	double sigma_,                                     ///< Écart-type initial ou calculé du noyau Gaussien.
	       theta_;                                     ///< Angle d’orientation (en radians) de la gaussienne si applicable.
	std::vector<double> pMin_,                         ///< Bornes inférieures autorisées pour les paramètres du fit.
	                    p_,                            ///< Valeurs courantes des paramètres du modèle.
	                    pMax_;                         ///< Bornes supérieures autorisées pour les paramètres du fit.
	static constexpr size_t N_ITER = 50;               ///< Nombre maximal d’itérations autorisées dans l’algorithme d’optimisation.
	static constexpr double EPSILON_1 = 1e-6,          ///< Seuil de convergence pour le changement de la fonction objectif.
	                        EPSILON_2 = 1e-6,          ///< Seuil de convergence pour le gradient.
	                        EPSILON_3 = 1e-5,          ///< Seuil de convergence pour la variation des paramètres.
	                        EPSILON_4 = 1e-2;          ///< Tolérance utilisée pour juger la stagnation ou instabilité.
	static constexpr double LAMBDA_0 = 1e-2,           ///< Valeur initiale du paramètre lambda dans l’algorithme de Levenberg-Marquardt.
	                        LAMBDA_EPSILON_MAX = 1e7,  ///< Valeur maximale autorisée pour lambda.
	                        LAMBDA_EPSILON_MIN = 1e-7, ///< Valeur minimale autorisée pour lambda.
	                        LAMBDA_UP = 11,            ///< Facteur multiplicatif lorsque la mise à jour est inefficace (augmentation de lambda).
	                        LAMBDA_DN = 9;             ///< Facteur de réduction de lambda lorsque la mise à jour est réussie.

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
	/// <param name="x">Coordonnée X estimée du centre initial.</param>
	/// <param name="y">Coordonnée Y estimée du centre initial.</param>
	/// <param name="sigmaX">Valeur initiale de la largeur gaussienne en X.</param>
	/// <param name="sigmaY">Valeur initiale de la largeur gaussienne en Y.</param>
	/// <param name="theta">Valeur initiale de la rotation (en radians, utilisée uniquement dans le mode 3D avec angle).</param>
	/// <param name="minI">Intensité minimale dans la ROI (servira à initialiser l'offset).</param>
	/// <param name="maxI">Intensité maximale dans la ROI (servira à initialiser l'amplitude).</param>
	void initP(double x, double y, double sigmaX, double sigmaY, double theta, double minI, double maxI);

	/// <summary>Calcule le coefficient de corrélation au carré (R²) entre deux vecteurs normalisés.
	/// Cela donne une mesure de la qualité d'ajustement entre les données observées et le modèle.</summary>
	/// <param name="dat">Vecteur des données observées (modifié en place).</param>
	/// <param name="hat">Vecteur des données ajustées (modifié en place).</param>
	/// <returns>Valeur de R², dans [0, 1].</returns>
	double computeRSS(std::vector<double>& dat, std::vector<double>& hat) const;

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
	/// <param name="alpha">Matrice approximative de Hessienne</param>
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
