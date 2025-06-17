///-------------------------------------------------------------------------------------------------
///
/// \file Utils.hpp
/// \brief Fonctions utiles pour la manipulation des matrices lors des ajustements.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <vector>

namespace palm_tracer {
constexpr size_t SPLINE_PATCH_SIDE   = 4,  ///< Taille d'un côté du patch spline (nombre de points en X ou Y).
                 SPLINE_PATCH_SQUARE = 16, ///< Nombre total de points dans un plan 2D du patch spline : SPLINE_PATCH_SIDE * SPLINE_PATCH_SIDE.
                 SPLINE_PATCH_SIZE   = 64; ///< Nombre total de coefficients dans un patch spline 3D : SPLINE_PATCH_SQUARE * SPLINE_PATCH_SIDE.

/// <summary> Structure représentant un pixel avec ses coordonnées (x, y) et une valeur scalaire associée (v).</summary>
struct Pixel {
	size_t x, ///< Coordonnée entière horizontale du pixel (axe X).
	       y; ///< Coordonnée entière verticale du pixel (axe Y).
	double v; ///< Valeur associée au pixel (intensité, poids, erreur, etc.).
};

//==================== PSF ====================
/// <summary>PSF (Point Spread Function) gaussienne isotrope (2D) centrée en (0,0), à sigma fixe.</summary>
/// <param name="offset">Valeur de fond (bruit ou signal de fond).</param>
/// <param name="i0">Intensité maximale de la gaussienne.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigma">Écart-type commun (isotrope).</param>
/// <returns>Valeur de la PSF en (x, y).</returns>
double PSF(double offset, double i0, double x, double y, double sigma);

/// <summary>PSF (Point Spread Function) gaussienne anisotrope 2D avec rotation (theta).</summary>
/// <param name="offset">Valeur de fond (bruit ou signal de fond).</param>
/// <param name="i0">Intensité maximale de la gaussienne.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigmaX">Écart-type selon l'axe X.</param>
/// <param name="sigmaY">Écart-type selon l'axe Y.</param>
/// <param name="theta">Angle de rotation (en radians).</param>
/// <returns>Valeur de la PSF en (x, y).</returns>
double PSF(double offset, double i0, double x, double y, double sigmaX, double sigmaY, double theta);

//==================== Matrix ====================
/// <summary>Dérivées partielles de la PSF 2D isotrope (centrée) par rapport aux paramètres.</summary>
/// <param name="out">Tableau de sortie (au moins 5 éléments).</param>
/// <param name="i0">Intensité maximale.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigma">Écart-type commun.</param>
/// <param name="has_sigma">True si sigma est un paramètre libre, false sinon.</param>
void Deriv(double* out, double i0, double x, double y, double sigma, bool has_sigma);

/// <summary>Dérivées partielles de la PSF 2D anisotrope avec rotation par rapport aux paramètres.</summary>
/// <param name="out">Tableau de sortie (au moins 7 éléments).</param>
/// <param name="i0">Intensité maximale.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigmaX">Écart-type selon X.</param>
/// <param name="sigmaY">Écart-type selon Y.</param>
/// <param name="theta">Angle de rotation.</param>
/// <param name="has_theta">True si theta est un paramètre libre, false sinon.</param>
void Deriv(double* out, double i0, double x, double y, double sigmaX, double sigmaY, double theta, bool has_theta);

/// <summary>Inversion d'une matrice carrée par la méthode de Gauss-Jordan.</summary>
/// <param name="in">Matrice à inverser (format vectorisé ligne par ligne, modifiée en place).</param>
/// <param name="out">Matrice résultat de l'inversion (même format).</param>
/// <param name="n"> Taille de la matrice (nxn). </param>
void InverseGaussJordan(std::vector<double>& in, std::vector<double>& out, size_t n);


/// <summary>Inversion de matrice (méthode LU (sans pivot) + systèmes triangulaires).</summary>
/// <param name="in">Matrice à inverser (format vectorisé ligne par ligne, modifiée en place).</param>
/// <param name="out">Matrice résultat de l'inversion (même format).</param>
/// <param name="n"> Taille de la matrice (nxn). </param>
void InverseLU(std::vector<double>& in, std::vector<double>& out, size_t n);


/// <summary>Calcule le centre de masse d'une image 2D carrée.</summary>
/// <param name="data">Données des pixels (valeurs d’intensité).</param>
/// <param name="n">Largeur/hauteur de l'image carrée (n x n).</param>
/// <param name="x">Coordonnées x du centre de masse.</param>
/// <param name="y">Coordonnées y du centre de masse.</param>
void CenterOfMass(const std::vector<Pixel>& data, size_t n, double& x, double& y);

/// <summary>Calcule les contributions d'une base de fonctions splines 3D (valeurs et dérivées).
/// Cette fonction calcule la valeur et les dérivées par rapport à x, y et z d'une fonction spline 3D basée sur des monômes (cx, cy, cz) à chaque point de l'espace (x, y, z).
/// Le résultat est stocké dans trois tableaux (valeur, dérivées partielles) de taille SPLINE_PATCH_SIDE^3.
/// </summary>
/// <param name="x">Coordonnée x (locale, typiquement entre 0 et 1).</param>
/// <param name="y">Coordonnée y (locale).</param>
/// <param name="z">Coordonnée z (locale).</param>
/// <param name="f">Référence vers le vecteur contenant les valeurs des fonctions de base.</param>
/// <param name="dxf">Référence vers le vecteur contenant les dérivées de f par rapport à x (∂f/∂x).</param>
/// <param name="dyf">Référence vers le vecteur contenant les dérivées de f par rapport à y (∂f/∂y).</param>
/// <param name="dzf">Référence vers le vecteur contenant les dérivées de f par rapport à z (∂f/∂z).</param>
void Delta3D(double x, double y, double z, std::vector<double>& f, std::vector<double>& dxf, std::vector<double>& dyf, std::vector<double>& dzf);

/// <summary>Réalise une décomposition de Cholesky d'une matrice symétrique définie positive.
///
/// Cette fonction effectue la décomposition d'une matrice réelle symétrique définie positive `a` de taille n×n selon la formule :
/// A = L * Lᵀ = Uᵀ * U.
///
/// Ici, L est triangulaire inférieure, U est triangulaire supérieure, et L = Uᵀ.
/// </summary>
/// <param name="a">Vecteur représentant la matrice A à décomposer (stockée en ligne principale).</param>
/// <param name="n">Dimension de la matrice A (A est de taille n x n).</param>
/// <param name="l">Référence vers le vecteur de sortie contenant la matrice L (triangulaire inférieure).</param>
/// <param name="u">Référence vers le vecteur de sortie contenant la matrice U (triangulaire supérieure).</param>
/// <returns><c>True</c> si la décomposition a réussi (matrice définie positive), <c>False</c> sinon.</returns>
/// <remarks>En cas d’échec (diagonale négative ou nulle), les vecteurs de sortie ne sont pas fiables.
/// L’utilisateur doit s’assurer que `l` et `u` sont préalloués de taille n x n.</remarks>
bool Cholesky(const std::vector<double>& a, size_t n, std::vector<double>& l, std::vector<double>& u);
}
