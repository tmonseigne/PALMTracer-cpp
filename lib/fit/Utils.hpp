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
constexpr size_t SPLINE_PATCH_SIDE   = 4,
                 SPLINE_PATCH_SQUARE = 16,
                 SPLINE_PATCH_SIZE   = 64;

//==================== PSF ====================
/// <summary>PSF gaussienne isotrope (2D) centrée en (0,0), à sigma fixe.</summary>
/// <param name="offset">Valeur de fond (background constant).</param>
/// <param name="i0">Amplitude de la gaussienne.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigma">Écart-type commun (isotrope).</param>
/// <returns>Valeur de la PSF en (x, y).</returns>
double PSF(double offset, double i0, double x, double y, double sigma);

/// <summary>PSF gaussienne anisotrope 2D avec rotation (theta).</summary>
/// <param name="offset">Valeur de fond.</param>
/// <param name="i0">Amplitude de la gaussienne.</param>
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
/// <param name="i0">Amplitude.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigma">Écart-type commun.</param>
/// <param name="has_sigma">True si sigma est un paramètre libre (ajouter dérivée), false sinon.</param>
void Deriv(double* out, double i0, double x, double y, double sigma, bool has_sigma);

/// <summary>Dérivées partielles de la PSF 2D anisotrope avec rotation par rapport aux paramètres.</summary>
/// <param name="out">Tableau de sortie (au moins 7 éléments).</param>
/// <param name="i0">Amplitude.</param>
/// <param name="x">Coordonnée x relative au centre.</param>
/// <param name="y">Coordonnée y relative au centre.</param>
/// <param name="sigmaX">Écart-type selon X.</param>
/// <param name="sigmaY">Écart-type selon Y.</param>
/// <param name="theta">Angle de rotation.</param>
/// <param name="has_theta">True si theta est un paramètre libre, false sinon.</param>
void Deriv(double* out, double i0, double x, double y, double sigmaX, double sigmaY, double theta, bool has_theta);

/// <summary>Inversion de matrice (méthode spécifique à l'implémentation, incomplète ici).</summary>
/// <param name="in">Matrice à inverser (format vectorisé ligne par ligne).</param>
/// <param name="out">Matrice résultat de l'inversion (même format).</param>
/// <param name="n"> Taille de la matrice (nxn). </param>
void Inverse(std::vector<double>& in, std::vector<double>& out, size_t n);

/// <summary></summary>
/// <param name="data"></param>
/// <param name="n"></param>
/// <param name="x"></param>
/// <param name="y"></param>
void CenterOfMass(const std::vector<double>& data, size_t n, double& x, double& y);

void Delta3D(double x, double y, double z, std::vector<double>& f, std::vector<double>& dxf, std::vector<double>& dyf, std::vector<double>& dzf);

bool Cholesky(const std::vector<double>& a, size_t n, std::vector<double>& l, std::vector<double>& u);
}
