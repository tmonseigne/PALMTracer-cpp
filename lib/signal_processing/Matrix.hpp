///-------------------------------------------------------------------------------------------------
///
/// \file Matrix.hpp
/// \brief Fonctions utiles pour la manipulation des matrices.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <algorithm>
#include <cmath>
#include <execution>
#include <numeric>
#include <ranges>
#include <vector>
#include "Structs.hpp"

namespace palm_tracer::matrix {
//=============================================
//==================== I/O ====================
//=============================================

/// <summary>Affiche le contenu d’un volume 3D (stocké en 1D) dans la console, plan par plan. </summary>
/// <param name="data">Données 1D représentant un volume 3D.</param>
/// <param name="size">Dimensions (X, Y, Z) du volume.</param>
void Print(const std::vector<double>& data, const Coord3D& size);

//========================================================
//==================== Basic Computes ====================
//========================================================

/// <summary>Vérifie si un vecteur contient au moins un NaN. </summary>
/// <param name="data">Vecteur à vérifier.</param>
/// <returns>True si au moins un élément est NaN, false sinon.</returns>
inline bool HasNaN(const std::vector<double>& data) { return std::ranges::any_of(data, [](const double x) { return std::isnan(x); }); }

/// <summary>Vérifie si un vecteur contient au moins une valeur non finie (NaN, infini). </summary>
/// <param name="data">Vecteur à vérifier.</param>
/// <returns>True si au moins un élément est NaN, false sinon.</returns>
inline bool IsValid(const std::vector<double>& data) { return std::ranges::any_of(data, [](const double x) { return !std::isfinite(x); }); }

/// <summary>Génère un masque booléen indiquant les éléments valides (non-NaN) du vecteur. </summary>
/// <param name="data">Vecteur d'entrée à tester.</param>
/// <returns>Vecteur booléen de même taille que data, avec true pour les valeurs non-NaN.</returns>
std::vector<bool> GetValidIndex(const std::vector<double>& data);

/// <summary>Compte le nombre d'éléments valides (ni NaN ni infini) dans un vecteur. </summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <returns>Nombre d'éléments valides.</returns>
size_t CountValid(const std::vector<double>& data);

/// <summary>Récupère l'indice du maximum du tableau.</summary>
/// <param name="data">Tableau.</param>
/// <returns>Indice du maximum.</returns>
inline size_t MaxId(const std::vector<double>& data) { return std::distance(data.begin(), std::ranges::max_element(data)); }

/// <summary>Valeur maximale dans un vecteur, en ignorant les NaN. </summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <returns>Valeur maximale non-NaN.</returns>
double Max(const std::vector<double>& data);

/// <summary>Somme des éléments d'un vecteur. </summary>
/// <param name="data">Vecteur à sommer.</param>
/// <returns>Somme des valeurs.</returns>
inline double Sum(const std::vector<double>& data) { return std::accumulate(data.begin(), data.end(), 0.0); }

/// <summary>Somme des valeurs valides (non-NaN, non-infini) d'un vecteur. </summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <returns>Somme des éléments valides.</returns>
double SumValid(const std::vector<double>& data);

/// <summary>Moyenne des éléments d’un vecteur. </summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <returns>Moyenne des éléments.</returns>
inline double Mean(const std::vector<double>& data) { return Sum(data) / static_cast<double>(data.size()); }

/// <summary>Moyenne des éléments valides (non-NaN, non-infini) d’un vecteur. </summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <returns>Moyenne des éléments valides, ou NaN si aucun n’est valide.</returns>
double MeanValid(const std::vector<double>& data);

/// <summary>Moyenne des éléments d'un vecteur selon un masque booléen. Les NaN sont ignorés automatiquement. </summary>
/// <param name="data">Données à moyenner.</param>
/// <param name="mask">Masque d'inclusion (doit avoir au moins la même taille que data).</param>
/// <returns>Moyenne des valeurs valides et masquées, ou NaN si aucune valeur valide.</returns>
double MeanWithMask(const std::vector<double>& data, const std::vector<bool>& mask);

/// <summary> Renvoie la médiane d’un vecteur de valeurs (en ignorant les NaN). </summary>
double Median(const std::vector<double>& data);

/// <summary> Renvoie un quantile donné (e.g., 0.75 pour le 75e percentile), en ignorant les NaN.</summary>
double Quantile(const std::vector<double>& data, double q);

/// <summary>Multiplie tous les éléments d'un vecteur par un facteur donné (modifie en place). </summary>
/// <param name="data">Vecteur à modifier.</param>
/// <param name="factor">Facteur multiplicatif.</param>
inline void Scale(std::vector<double>& data, const double factor) { for (auto& x : data) { x *= factor; } }

/// <summary>Produit scalaire entre deux vecteurs. </summary>
/// <param name="a">Premier vecteur.</param>
/// <param name="b">Deuxième vecteur.</param>
/// <returns>Somme des produits élément par élément.</returns>
inline double DotProduct(const std::vector<double>& a, const std::vector<double>& b) { return std::transform_reduce(a.begin(), a.end(), b.begin(), 0.0); }

/// <summary>Calcule la différence quadratique moyenne entre deux vecteurs, en ignorant les NaN. </summary>
/// <param name="a">Premier vecteur.</param>
/// <param name="b">Deuxième vecteur (même taille).</param>
/// <returns>La différence quadratique moyenne.</returns>
double RmsDifference(const std::vector<double>& a, const std::vector<double>& b);

/// <summary>Décalage vers 0 et normalisation de la matrice.
/// Si `shift` est vrai, la fonction soustrait le minimum à tous les éléments, puis divise par le maximum pour obtenir un tableau entre 0 et 1.
/// Si `shift` est faux, seule la division par le maximum est appliquée. </summary>
/// <param name="a">Matrice à normaliser (modifié en place).</param>
/// <param name="shift">Autorise le décallage vers 0 du minimum.</param>
void Normalisation(std::vector<double>& a, bool shift = true);

/// <summary>Coupe les bords du volume sur les dimensions spécifiées.
///
/// Cette fonction supprime une couche de voxels sur les bords de chaque dimension (x, y, z) spécifiée à true.
/// Elle modifie le vecteur de données en place et met à jour les dimensions.
/// </summary>
/// <param name="data">Données 1D représentant un volume 3D (modifiées en place).</param>
/// <param name="size">Dimensions initiales du volume (modifiées en place).</param>
/// <param name="cropX">Coupe la première et dernière colonne (x).</param>
/// <param name="cropY">Coupe la première et dernière ligne (y).</param>
/// <param name="cropZ">Coupe le premier et dernier plan (z).</param>
void Crop(std::vector<double>& data, Coord3D& size, bool cropX = true, bool cropY = true, bool cropZ = true);

/// <summary> Extrait un sous-volume 3D depuis un volume d'entrée 1D donné des indices spécifiques. </summary>
/// <param name="data">Données 1D représentant un volume 3D.</param>
/// <param name="size">Dimensions initiales du volume.</param>
/// <param name="xrange">Indices en X.</param>
/// <param name="yrange">Indices en Y.</param>
/// <param name="zrange">Indices en Z.</param>
/// <returns>Volume extrait (1D, de taille xr.size() * yr.size() * zr.size()).</returns>
std::vector<double> ExtractSubVolume(const std::vector<double>& data, const Coord3D& size, const std::vector<size_t>& xrange,
                                     const std::vector<size_t>& yrange, const std::vector<size_t>& zrange);

//======================================================
//==================== Manipulation ====================
//======================================================

/// <summary>Calcule le centre de masse d'une image 2D carrée.</summary>
/// <param name="data">Données des pixels (valeurs d’intensité).</param>
/// <param name="n">Largeur/hauteur de l'image carrée (n x n).</param>
/// <param name="x">Coordonnées x du centre de masse.</param>
/// <param name="y">Coordonnées y du centre de masse.</param>
void CenterOfMass(const std::vector<Pixel>& data, size_t n, double& x, double& y);

/// <summary>Réalise une décomposition de Cholesky d'une matrice symétrique définie positive.
///
/// Cette fonction effectue la décomposition d'une matrice réelle symétrique définie positive `A` de taille nxn selon la formule :
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
} // namespace palm_tracer::matrix
