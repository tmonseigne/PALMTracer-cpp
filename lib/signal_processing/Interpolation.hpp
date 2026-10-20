///-------------------------------------------------------------------------------------------------
///
/// \file Interpolation.hpp
/// <summary>Outils d'interpolation 3D : zoom, décalage et estimation de shift subpixel.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 02/07/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// <remarks>Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.</remarks>
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include "Structs.hpp"


namespace palm_tracer {
/// <summary>Méthode d'interpolation.</summary>
enum class interpolation_mode : uint8_t {
	trilinear,    ///< Interpolation trilineaire (rapide).
	cubic,        ///< Interpolation Cubique classique (Catmull-Rom).
	cubic_weight, ///< Interpolation spline cubique pondérée (a = -0.5), proche de MATLAB.
	spline        ///< Interpolation spline spline B-cubique (matlab-like, spline contigu).
};

/// <summary>Classe gérant l'interpolation 3D (zoom, shift, alignement subpixel).</summary>
class Interpolation3D {
public:
	interpolation_mode mode_; ///< Mode d'interpolation à utiliser.
	bool cropBound_;          ///< Si activé : toute coordonnée hors volume renvoie 0.0 (sinon accès réflexif).

	/// <summary>Constructeur de la classe <c>Interpolation3D</c>.</summary>
	/// <param name="mode">Mode d'interpolation à utiliser.</param>
	/// <param name="cropBound">Active la suppression hors-bord.</param>
	explicit Interpolation3D(const interpolation_mode mode = interpolation_mode::cubic_weight, const bool cropBound = false)
		: mode_(mode), cropBound_(cropBound) { }

	/// <summary>Applique un zoom isotrope sur un volume 3D avec interpolation selon le mode.</summary>
	/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
	/// <param name="level">Facteur de zoom entier (> 0). Un niveau 1 insère un point entre chaque valeur initiale, un niveau 2 en insère 3, etc.
	/// Soit \f$2^l - 1\f$ points interpolés entre chaque voxel.</param>
	/// <param name="newDim">Taille du volume de sortie, calculée comme :
	/// - `(size.x - 1) * 2^level + 1`
	/// - `(size.y - 1) * 2^level + 1`
	/// - `(size.z - 1) * 2^level + 1`
	/// Ce comportement imite fidèlement le comportement de MATLAB `interp3(M, level)`.</param>
	/// <returns>Volume agrandi sous forme linéarisée avec interpolation.</returns>
	[[nodiscard]] std::vector<double> Zoom(const std::vector<double>& volume, const Vec3I& dim, size_t level, Vec3I& newDim) const;

	/// <summary>Applique un décalage (shift) flottant sur un volume 3D avec interpolation.</summary>
	/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
	/// <param name="shift">Décalage flottant (x, y, z).</param>
	/// <returns>Volume décalé, même taille que l'entrée.</returns>
	[[nodiscard]] std::vector<double> Shift(const std::vector<double>& volume, const Vec3I& dim, const Vec3D& shift) const;

	/// <summary>Estime le décalage 3D subpixel entre deux volumes par corrélation croisée en haute résolution.
	///
	/// Cette méthode aligne un volume cible `target` sur un volume de référence `ref` en utilisant :
	/// - un suréchantillonnage du volume (par interpolation trilineaire ou cubique),
	/// - une transformée de Fourier 3D (FFT),
	/// - une corrélation croisée dans le domaine fréquentiel (`ref * conj(target)`),
	/// - un recentrage (`fftshift`) pour que le maximum apparaisse au centre,
	/// - et une localisation du pic de corrélation par interpolation fine (`SubpixelMaxInterp3D`).
	///
	/// Le décalage estimé est retourné avec une précision subpixel.</summary>
	/// <param name="ref">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="target">Volume cible à aligner (entrée), 3D en mémoire row-major (x, y, z).</param>
	/// <param name="dim">Dimensions des volumes sources (identiques pour les deux).</param>
	/// <returns>Décallage 3D subpixel estimé.</returns>
	[[nodiscard]] Vec3D GetShift(const std::vector<double>& ref, const std::vector<double>& target, const Vec3I& dim) const;

	/// <summary>Recherche la position subvoxel du maximum local d'un volume 3D par interpolation locale.</summary>
	/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="dim">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
	/// <param name="pos">Position entière estimée du maximum [x, y, z].</param>
	/// <param name="delta">Pas d'échantillonnage dans l'interpolation locale.</param>
	/// <param name="w">Demi-taille de la fenêtre (1 ⇒ cube 3x3x3).</param>
	/// <returns>Vecteur (x, y, z) de décalage local subvoxel.</returns>
	static Vec3D SubpixelMax(const std::vector<double>& volume, const Vec3I& dim, const Vec3D& pos, double delta = 0.05, double w = 1.0);

	/// <summary>Réalise un suréchantillonnage isotrope par un facteur flottant.
	/// Équivalent à interp3(..., 'method') avec Xq = 1:1/factor:N.
	/// </summary>
	/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="dim">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
	/// <param name="factor">Facteur d'agrandissement (> 1.0), ex: 3.0 pour 3x.</param>
	/// <param name="newDim">Dimensions du volume suréchantillonné en sortie.</param>
	/// <returns>Volume interpolé à la taille souhaitée.</returns>
	[[nodiscard]] std::vector<double> Upsample(const std::vector<double>& volume, const Vec3I& dim, double factor, Vec3I& newDim) const;

private:
	/// <summary>Teste si une coordonnée réelle est hors des bornes d'une dimension.</summary>
	/// <param name="value">Coordonnée.</param>
	/// <param name="dim">Dimensions du volume (width=X, height=Y, depth=Z).</param>
	/// <returns><c>True</c> Si la valeur est hors cadre, <c>False</c> sinon.</returns>
	static bool isOut(const double value, const size_t dim) { return value < 0 || value >= static_cast<double>(dim); }

	/// <summary>Teste si un point 3D est hors volume.</summary>
	/// <param name="pos">Position du point.</param>
	/// <param name="dim">Dimensions du volume (width=X, height=Y, depth=Z).</param>
	/// <returns><c>True</c> Si la position est hors cadre, <c>False</c> sinon.</returns>
	static bool isOut(const Vec3D& pos, const Vec3I& dim) { return isOut(pos.x, dim.x) || isOut(pos.y, dim.y) || isOut(pos.z, dim.z); }

	/// <summary>Effectue une interpolation cubique 1D sur quatre valeurs. Basé sur la méthode de Catmull-Rom (forme cardinale).</summary>
	/// <param name="v">Tableau de 4 valeurs consécutives.</param>
	/// <param name="x">Position relative entre les points (dans [0, 1]).</param>
	/// <returns>Valeur interpolée à la position x.</returns>
	static double cubic1D(const std::array<double, 4>& v, double x);

	/// <summary>Effectue une interpolation cubique 1D (formule pondérée \f$a = -0.5\f$).</summary>
	/// <param name="x">Position relative entre les points (dans [0, 1]).</param>
	/// <param name="a">Coefficient de la spline Cubique.</param>
	/// <returns>Valeur interpolée à la position x.</returns>
	static double cubicWeight1D(double x, double a = -0.5);

	/// <summary>Interpolation cubique trivariée dans un volume 3D (64 voxels pris en compte, 4x4x4).</summary>
	/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="dim">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
	/// <param name="pos">Coordonnées flottantes où interpoler.</param>
	/// <returns>Valeur interpolée au point donné, ou 0.0 si hors volume.</returns>
	static double cubic3D(const std::vector<double>& volume, const Vec3I& dim, const Vec3D& pos);

	/// <summary>Interpolation cubique trivariée dans un volume 3D (64 voxels pris en compte, 4x4x4).</summary>
	/// <param name="volume">Vecteur linéarisé représentant un volume 3D en mémoire ligne-colonne-plan (X-Y-Z).</param>
	/// <param name="dim">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
	/// <param name="pos">Coordonnées flottantes où interpoler.</param>
	/// <returns>Valeur interpolée au point donné, ou 0.0 si hors volume.</returns>
	static double cubicWeight3D(const std::vector<double>& volume, const Vec3I& dim, const Vec3D& pos);

	/// <summary>Effectue une interpolation trilineaire dans un volume 3D (rapide, linéaire dans les 3 axes).</summary>
	/// <param name="volume">Vecteur linéarisé représentant le volume 3D (X-Y-Z) en mémoire row-major.</param>
	/// <param name="dim">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
	/// <param name="pos">Position flottante (x, y, z) où calculer l'interpolation.</param>
	/// <returns>Valeur interpolée au point donné, ou 0.0 si hors volume.</returns>
	static double trilinear(const std::vector<double>& volume, const Vec3I& dim, const Vec3D& pos);

	/// <summary>Interpolation spline B-cubique avec coefficients spline préalablement calculés.</summary>
	/// <param name="coeffs">Coefficients de la spline pour le volume.</param>
	/// <param name="dim">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
	/// <param name="pos">Position flottante (x, y, z) où calculer l'interpolation.</param>
	/// <returns>Valeur interpolée au point donné, ou 0.0 si hors volume.</returns>
	static double spline3D(const std::vector<double>& coeffs, const Vec3I& dim, const Vec3D& pos);
};
} // namespace palm_tracer
