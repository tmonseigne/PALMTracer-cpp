///-------------------------------------------------------------------------------------------------
///
/// \file PSF.hpp
/// \brief Fonctions utiles pour le calcul des PSF.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 01/07/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <vector>

namespace palm_tracer {
/// <summary>PSF (Point Spread Function) Gaussienne **isotrope** en 2D (centre implicite à l'origine).\n
/// Modèle : \f$\mathrm{PSF}(x,y)=\text{offset}+I_0\exp\!\left(-\dfrac{x^2+y^2}{2\sigma^2}\right)\f$.</summary>
/// <param name="offset">Niveau de fond ajouté (même unité que l’intensité).</param>
/// <param name="i0">Amplitude crête \f$I_0\f$ (intensité du pic).</param>
/// <param name="x">Abscisse **relative** (en pixels) : \f$x=x_\mathrm{obs}-x_0\f$.</param>
/// <param name="y">Ordonnée **relative** (en pixels) : \f$y=y_\mathrm{obs}-y_0\f$.</param>
/// <param name="sigma">Écart-type isotrope \f$\sigma>0\f$ (en pixels).</param>
/// <returns>Valeur scalaire de la PSF au point \f$(x, y)\f$.</returns>
/// <remarks>
/// - Les coordonnées \f$(x, y)\f$ doivent être **centrées** par rapport au paramètre de position \f$(x_0, y_0)\f$ du modèle.
/// - Stabilité numérique : éviter \f$\sigma \to 0\f$ (risque de sous-/sur-flux via l’exponentielle).
/// </remarks>
double PSF(double offset, double i0, double x, double y, double sigma);


/// <summary>PSF (<i>Point Spread Function</i>) Gaussienne **anisotrope et orientée** en 2D (centre implicite à l'origine) avec rotation (theta).
/// Modèle : rotation de \f$\theta\f$ puis dilatation anisotrope \f$(\sigma_x, \sigma_y)\f$
/// \f$(x',y')=R(-\theta, x, y),\quad \mathrm{PSF}(x,y)=\text{offset}+I_0\exp\!\big(-\,(x'/\sigma_x)^2-\,(y'/\sigma_y)^2\big)\f$.</summary>
/// <param name="offset">Niveau de fond (même unité que l’intensité).</param>
/// <param name="i0">Amplitude crête \f$I_0\f$.</param>
/// <param name="x">Abscisse **relative** (pixels) : \f$x=x_\mathrm{obs}-x_0\f$.</param>
/// <param name="y">Ordonnée **relative** (pixels) : \f$y=y_\mathrm{obs}-y_0\f$.</param>
/// <param name="sigmaX">\f$\sigma_x > 0\f$ (pixels).</param>
/// <param name="sigmaY">\f$\sigma_y > 0\f$ (pixels).</param>
/// <param name="theta">Angle d’orientation \f$\theta\f$ (radians, convention math. directe).</param>
/// <returns> Valeur scalaire de la PSF au point \f$(x,y)\f$. </returns>
/// <remarks>
/// - Convention de rotation : \f$R(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{pmatrix}\f$,
///   on applique \f$(x',y')=R(-\theta)\,(x,y)\f$ (équivalent aux formules implantées).
/// - Les coordonnées \f$(x,y)\f$ sont **centrées** (cf. <see cref="GaussianFit::psf"/>).
/// - Stabilité : éviter \f$\sigma_x, \sigma_y \to 0\f$.
/// </remarks>
double PSF(double offset, double i0, double x, double y, double sigmaX, double sigmaY, double theta);

/// <summary>Gradient de la PSF **isotrope** (modèle 2D circulaire) en coordonnées **centrées**.
/// Le vecteur de sortie est :
/// \f[
/// \left[\dfrac{\partial}{\partial\text{offset}},\,
///               \dfrac{\partial}{\partial I_0},\,
///               \dfrac{\partial}{\partial x_0},\,
///               \dfrac{\partial}{\partial y_0},\,
///               \dfrac{\partial}{\partial \sigma}\right]
/// \f]
/// (le dernier terme n’est écrit que si <paramref name="has_sigma"/> vaut <c>true</c>).
/// </summary>
/// <param name="out">Vecteur de sortie (taille ≥ 5 si <paramref name="has_sigma"/>).</param>
/// <param name="i0">Amplitude \f$I_0\f$.</param>
/// <param name="x">Abscisse **relative** (pixels) : \f$x=x_\mathrm{obs}-x_0\f$.</param>
/// <param name="y">Ordonnée **relative** (pixels) : \f$y=y_\mathrm{obs}-y_0\f$.</param>
/// <param name="sigma">\f$\sigma>0\f$ (pixels).</param>
/// <param name="has_sigma">Si <c>true</c>, calcule aussi \f$\partial\mathrm{PSF}/\partial \sigma\f$.</param>
/// <remarks>
/// - Rappels analytiques (avec \f$E=\exp\big(-\dfrac{x^2+y^2}{2\sigma^2}\big)\f$ et \f$P=I_0 E\f$) : \n
///   \f[
///   \dfrac{\partial}{\partial\text{offset}}=1,\quad
///   \dfrac{\partial}{\partial I_0} = E,\quad
///   \dfrac{\partial}{\partial x_0} = P\,\dfrac{x}{\sigma^2},\quad
///   \dfrac{\partial}{\partial y_0} = P\,\dfrac{y}{\sigma^2},\quad
///   \dfrac{\partial}{\partial \sigma} =  P\,\dfrac{x^2+y^2}{\sigma^3}
///   \f]
/// - Les signes apparents liés à \f$\dfrac{\partial}{\partial x_0}\f$ et \f$\dfrac{\partial}{\partial y_0}\f$ sont déjà pris en compte
///   via l’utilisation des coordonnées **centrées** \f$(x,y)=(x_\mathrm{obs}-x_0,y_\mathrm{obs}-y_0)\f$.
/// </remarks>
void Deriv(std::vector<double>& out, double i0, double x, double y, double sigma, bool has_sigma);

/// <summary>Gradient de la PSF **anisotrope orientée** en 2D, coordonnées **centrées**.
/// Le vecteur de sortie est :
/// \f[
/// \left[\dfrac{\partial}{\partial\text{offset}}, \,
///             \dfrac{\partial}{\partial I_0}, \,
///             \dfrac{\partial}{\partial x_0}, \,
///             \dfrac{\partial}{\partial y_0}, \,
///             \dfrac{\partial}{\partial \sigma_x}, \,
///             \dfrac{\partial}{\partial \sigma_y}, \,
///             \dfrac{\partial}{\partial \theta}\right]
/// \f]
/// (le dernier terme n’est écrit que si <paramref name="has_theta"/> vaut <c>true</c>).
/// </summary>
/// <param name="out">Vecteur de sortie (taille ≥ 7 si <paramref name="has_theta"/>).</param>
/// <param name="i0">Amplitude \f$I_0\f$.</param>
/// <param name="x">Abscisse **relative** (pixels) : \f$x=x_\mathrm{obs}-x_0\f$.</param>
/// <param name="y">Ordonnée **relative** (pixels) : \f$y=y_\mathrm{obs}-y_0\f$.</param>
/// <param name="sigmaX">\f$\sigma_x>0\f$ (pixels).</param>
/// <param name="sigmaY">\f$\sigma_y>0\f$ (pixels).</param>
/// <param name="theta">Angle \f$\theta\f$ (radians).</param>
/// <param name="has_theta">Si <c>true</c>, calcule aussi \f$\partial/\partial \theta\f$.</param>
/// <remarks>
/// - On utilise la même convention que <see cref="PSF(double, double, double, double, double, double, double)"/> :
///   \f$(x', y')=R(-\theta)\, (x, y)\f$ puis \f$\exp(-((x'/\sigma_x)^2+(y'/\sigma_y)^2))\f$.
/// - Les expressions codées correspondent aux dérivées analytiques (y compris termes croisés liés à \f$\theta\f$).
/// - Pour la stabilité, éviter \f$\sigma_x,\sigma_y\to 0\f$.
/// </remarks>
void Deriv(std::vector<double>& out, double i0, double x, double y, double sigmaX, double sigmaY, double theta, bool has_theta);
} // namespace palm_tracer
