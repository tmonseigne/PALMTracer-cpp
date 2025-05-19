///-------------------------------------------------------------------------------------------------
///
/// \file ATrous.hpp
/// \brief  Classe pour la décomposition en ondelettes 2D selon l’algorithme "à trous".
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <vector>

/// <summary> Classe pour la décomposition en ondelettes 2D selon l’algorithme "à trous".</summary>
/// <remarks>
/// Cette classe réalise la transformation en ondelettes sur une image 2D, notamment avec une décomposition multi-échelle utilisant des filtrages spline B3.
/// Elle gère la mémoire interne, les calculs de convolutions et la gestion des plans.
/// </remarks>
class ATrous
{
public:
	size_t height_,           ///< Hauteur de l'image analysée.
	       width_,            ///< Largeur de l'image analysée.
	       size_;             ///< Taille totale (hauteur × largeur).
	std::vector<double> tmp_; ///< Vecteur temporaire utilisé pour les calculs intermédiaires lors de la décomposition.
	std::vector<std::vector<double>> wDecomp_;
	///< Résultat de la décomposition en ondelettes : tableau des plans (chaque plan est un vecteur 2D de coefficients).

	/// <summary>Constructeur de la classe <c>ATrous</c>.</summary>
	ATrous() = default;

	/// <summary>Destructeur de la classe <c>ATrous</c>.</summary>
	~ATrous() = default;

	/// <summary>Applique la décomposition en ondelettes 2D par l’algorithme à trous sur une image d’entrée.</summary>
	/// <param name="image">Pointeur vers l’image 2D d’entrée (tableau d'unsigned short).</param>
	/// <param name="height">Hauteur de l’image.</param>
	/// <param name="width">Largeur de l’image.</param>
	/// <param name="planesToUse">Tableau booléen précisant quels plans doivent être calculés (doit être alloué).</param>
	/// <param name="nbPlane">Nombre de plans à calculer (sera ajusté si nécessaire).</param>
	/// <returns>Nombre de plans effectivement calculés, ou W_SIZE_ERR en cas d’erreur.</returns>
	int ATrousDecomposition2D(const unsigned short* image, size_t height, size_t width, bool* planesToUse, size_t nbPlane);

private:
	/// <summary>Calcule 2^val.</summary>
	/// <param name="val">Exposant (taille non signée)</param>
	/// <returns>Renvoie 2^val (2 élevé à la puissance val)</returns>
	/// <remarks>
	/// Cette fonction est constexpr, donc évaluée à la compilation si possible.
	/// Le mot-clé [[nodiscard]] indique que le résultat ne doit pas être ignoré.
	/// </remarks>
	[[nodiscard]] static constexpr size_t twoToThe(const size_t val) noexcept { return (val == 0) ? 1 : (static_cast<size_t>(1) << val); }

	/// <summary>Calcule l'offset dans un tableau 1D à partir d'une position 2D avec gestion des débordements par clamp.</summary>
	/// <param name="i">Indice de ligne (peut être négatif, sera borné entre 0 et height-1).</param>
	/// <param name="j">Indice de colonne (peut être négatif, sera borné entre 0 et width-1).</param>
	/// <returns>Offset 1D dans le tableau.</returns>
	[[nodiscard]] size_t offsetTest(int i, int j) const;

	/// <summary>Calcule l'offset dans un tableau 1D à partir d'une position 2D sans vérification.</summary>
	/// <param name="i">Indice de ligne (doit être valide).</param>
	/// <param name="j">Indice de colonne (doit être valide).</param>
	/// <returns>Offset 1D dans le tableau.</returns>
	[[nodiscard]] size_t offset(const size_t i, const size_t j) const { return i * width_ + j; }

	/// <summary>Vérifie si le nombre de plans spécifié est compatible avec les dimensions de l'image, en adaptant si nécessaire.</summary>
	/// <param name="planesToUse">Tableau booléen des plans à utiliser (doit être alloué avec au moins 9 cases).</param>
	/// <param name="nbPlane">Nombre de plans souhaité (sera décrémenté si trop grand).</param>
	/// <returns>Nombre réel de plans utilisables, ou W_SIZE_ERR en cas d’erreur.</returns>
	int testSizes2D(bool* planesToUse, size_t& nbPlane) const;

	/// <summary>Applique le filtre B3-spline 2D à l’image d’entrée pour générer une image de sortie à une certaine échelle.</summary>
	/// <param name="in">Image d’entrée (provenant d’un niveau de décomposition précédent).</param>
	/// <param name="out">Image de sortie filtrée.</param>
	/// <param name="plane">Indice du plan (échelle) utilisé pour ajuster les pas de la convolution.</param>
	void b3SplineFiltering2D(const std::vector<double>& in, std::vector<double>& out, size_t plane);
};
