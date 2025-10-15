///-------------------------------------------------------------------------------------------------
///
/// \file Matrix.hpp
/// \brief Fonctions utiles pour la manipulation des matrices et volumes 3D.
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

/// <summary>Affiche le contenu d’un volume 3D (stocké en 1D) dans la console, plan par plan.
/// Format lisible avec délimitation par lignes et colonnes.</summary>
/// <param name="data">Vecteur de données 1D (flattened) représentant un volume 3D.</param>
/// <param name="dim">Dimensions du volume (x, y, z).</param>
/// <remarks>Les valeurs sont affichées plan par plan avec une mise en forme tabulaire.</remarks>
void PrintVolume(const std::vector<double>& data, const Vec3I& dim);

/// <summary>Affiche le contenu d’une matrice 2D (stocké en 1D) dans la console.
/// Format lisible avec délimitation par lignes et colonnes.</summary>
/// <param name="data">Vecteur de données 1D (flattened) représentant une matrice 2D.</param>
/// <param name="row">Nombre de lignes de la matrice.</param>
/// <param name="col">Nombre de colonnes de la matrice.</param>
void PrintMatrix(const std::vector<double>& data, size_t row, size_t col);

/// <summary>Affiche le contenu d’une matrice 2D dans la console.
/// Format lisible avec délimitation par lignes et colonnes.</summary>
/// <param name="data">Matrice 2D.</param>
void PrintMatrix(const std::vector<std::vector<double>>& data);

/// <summary>Affiche le contenu d’un vecteur. </summary>
/// <param name="data">Vecteur de données 1D.</param>
void PrintVector(const std::vector<double>& data);

/// <summary>Sérialise <c>data</c> dans <c>output</c> avec le numéro de ligne, au format <c>[L, RowIndex, valeurs..., ..., -1]</c>.\n
/// L = 1 + row.size(). La fin du flux est marquée par -1.0 (L < 0 côté parsing).
/// </summary>
/// <param name="data">Tableau 2D irrégulier.</param>
/// <param name="output">Pointeur vers le buffer de sortie.</param>
/// <returns>Le nombre d'éléments écrits dans l'output (en comptant la valeur sentinelle de -1).</returns>
size_t SendIrregularArray(const std::vector<std::vector<double>>& data, double* output);

//========================================================
//==================== Basic Computes ====================
//========================================================

/// <summary>Vérifie si un vecteur contient au moins un NaN.</summary>
/// <param name="data">Vecteur à vérifier.</param>
/// <returns>True si au moins un élément est NaN, false sinon.</returns>
inline bool HasNaN(const std::vector<double>& data) { return std::ranges::any_of(data, [](const double x) { return std::isnan(x); }); }

/// <summary>Vérifie si un vecteur contient au moins une valeur non finie (NaN ou infini).</summary>
/// <param name="data">Vecteur à vérifier.</param>
/// <returns>True si au moins un élément est non fini (NaN ou infini), false sinon.</returns>
inline bool HasInvalid(const std::vector<double>& data) { return std::ranges::any_of(data, [](const double x) { return !std::isfinite(x) || std::isnan(x); }); }

/// <summary>Vérifie si un vecteur ne contient que des valeurs non finies (NaN ou infini).</summary>
/// <param name="data">Vecteur à vérifier.</param>
/// <returns>True si tous les éléments sont non finies (NaN ou infini), false sinon.</returns>
inline bool IsInvalid(const std::vector<double>& data) { return std::ranges::all_of(data, [](const double x) { return !std::isfinite(x) || std::isnan(x); }); }

/// <summary>Génère un masque booléen indiquant les éléments valides (finis et non-NaN) du vecteur.</summary>
/// <param name="data">Vecteur d'entrée à tester.</param>
/// <returns>Vecteur booléen de même taille que data, avec true pour les éléments valides.</returns>
std::vector<bool> GetValidIndex(const std::vector<double>& data);

/// <summary>Compte le nombre d'éléments valides (ni NaN, ni infini) dans un vecteur.</summary>
/// <param name="data">Vecteur d'entrée à tester.</param>
/// <returns>Nombre d'éléments valides.</returns>
size_t CountValid(const std::vector<double>& data);

/// <summary>Récupère l'indice du maximum du tableau (sans vérification de validité).</summary>
/// <param name="data">Vecteur d'entrée à tester.</param>
/// <returns>Indice de la valeur maximale.</returns>
/// <remarks>Les NaN ne sont pas ignorés.</remarks>
inline size_t MaxId(const std::vector<double>& data) { return std::distance(data.begin(), std::ranges::max_element(data)); }

/// <summary>Valeur maximale dans un vecteur, en ignorant les NaN et infinis.</summary>
/// <param name="data">Vecteur d'entrée à tester.</param>
/// <returns>Valeur maximale valide, ou -infini si aucune n'est valide.</returns>
/// <remarks>Les valeurs non finies sont ignorées.</remarks>
double Max(const std::vector<double>& data);

/// <summary>Somme des éléments d'un vecteur (sans vérification de validité).</summary>
/// <param name="data">Vecteur d'entrée à sommer.</param>
/// <returns>Somme des valeurs.</returns>
/// <remarks>Les NaN ne sont pas ignorés.</remarks>
inline double Sum(const std::vector<double>& data) { return std::accumulate(data.begin(), data.end(), 0.0); }

/// <summary>Somme des valeurs valides (ni NaN, ni infini) d'un vecteur.</summary>
/// <param name="data">Vecteur d'entrée à sommer.</param>
/// <returns>Somme des éléments valides.</returns>
double SumValid(const std::vector<double>& data);

/// <summary>Moyenne des éléments d’un vecteur (sans vérification de validité).</summary>
/// <param name="data">Vecteur d'entrée à moyenner.</param>
/// <returns>Moyenne des éléments.</returns>
inline double Mean(const std::vector<double>& data) { return Sum(data) / static_cast<double>(data.size()); }

/// <summary>Moyenne des éléments valides  (ni NaN, ni infini) d’un vecteur.</summary>
/// <param name="data">Vecteur d'entrée à moyenner.</param>
/// <returns>Moyenne des éléments valides, ou NaN si aucun n’est valide.</returns>
double MeanValid(const std::vector<double>& data);

/// <summary>Moyenne des éléments d'un vecteur selon un masque booléen. Les NaN sont ignorés automatiquement.</summary>
/// <param name="data">Vecteur d'entrée à moyenner.</param>
/// <param name="mask">Masque d'inclusion (doit avoir au moins la même taille que data).</param>
/// <returns>Moyenne des valeurs valides et masquées, ou NaN si aucune valeur valide.</returns>
double MeanWithMask(const std::vector<double>& data, const std::vector<bool>& mask);

/// <summary>Renvoie la médiane d’un vecteur, en ignorant les NaN et infinis.</summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <returns>Médiane, ou NaN si aucune valeur valide.</returns>
double Median(const std::vector<double>& data);

/// <summary>Calcule un quantile donné dans un vecteur, en ignorant les NaN et infinis.</summary>
/// <param name="data">Vecteur d'entrée.</param>
/// <param name="q">Quantile entre 0.0 et 1.0.</param>
/// <returns>Valeur du quantile, ou NaN si aucune valeur valide.</returns>
double Quantile(const std::vector<double>& data, double q);

/// <summary>Ajoute une constante à tous les éléments du vecteur.</summary>
/// <param name="data">Vecteur modifié en place.</param>
/// <param name="value">Constante à ajouter.</param>
inline void Add(std::vector<double>& data, const double value) { for (auto& x : data) { x += value; } }

/// <summary>Multiplie tous les éléments du vecteur par une constante.</summary>
/// <param name="data">Vecteur modifié en place.</param>
/// <param name="factor">Facteur multiplicatif.</param>
inline void Scale(std::vector<double>& data, const double factor) { for (auto& x : data) { x *= factor; } }

/// <summary>Produit scalaire classique entre deux vecteurs.</summary>
/// <param name="a">Premier vecteur.</param>
/// <param name="b">Deuxième vecteur.</param>
/// <returns>Somme des produits élément par élément.</returns>
inline double DotProduct(const std::vector<double>& a, const std::vector<double>& b) { return std::transform_reduce(a.begin(), a.end(), b.begin(), 0.0); }

/// <summary>Calcule la différence quadratique moyenne entre deux vecteurs, en ignorant les NaN.</summary>
/// <param name="a">Premier vecteur.</param>
/// <param name="b">Deuxième vecteur (même taille).</param>
/// <returns>RMS des différences valides, ou NaN si aucun élément valide.</returns>
double RmsDifference(const std::vector<double>& a, const std::vector<double>& b);

/// <summary>Normalise une matrice en la centrant puis en divisant par la somme totale.
/// Si `shift` est vrai, la fonction soustrait le minimum à tous les éléments et divise par le maximum pour obtenir un tableau entre 0 et 1.
/// Si `shift` est faux, seule la division par le maximum est appliquée.</summary>
/// <param name="a">Vecteur à normaliser (modifié en place).</param>
/// <param name="shift">Si true, retire le minimum avant normalisation.</param>
/// <remarks>La somme des éléments résultants est 1 (sauf cas dégénéré).</remarks>
void Normalisation(std::vector<double>& a, bool shift = true);

/// <summary>Supprime les bords du volume selon les axes sélectionnés.
///
/// Cette fonction supprime une couche de voxels sur les bords de chaque dimension (x, y, z) spécifiée à true.
/// Elle modifie le vecteur de données en place et met à jour les dimensions.
/// </summary>
/// <param name="data">Vecteur de données 1D (flattened) représentant un volume 3D.</param>
/// <param name="size">Dimensions initiales du volume (modifiées en place).</param>
/// <param name="cropX">Coupe la première et la dernière colonne (x).</param>
/// <param name="cropY">Coupe la première et la dernière ligne (y).</param>
/// <param name="cropZ">Coupe le premier et le dernier plan (z).</param>
void Crop(std::vector<double>& data, Vec3I& size, bool cropX = true, bool cropY = true, bool cropZ = true);

/// <summary>Extrait un sous-volume 3D donné des plages d'indices spécifiées.</summary>
/// <param name="data">Données 1D représentant un volume 3D.</param>
/// <param name="dim">Dimensions initiales du volume.</param>
/// <param name="xrange">Indices en X.</param>
/// <param name="yrange">Indices en Y.</param>
/// <param name="zrange">Indices en Z.</param>
/// <returns>Sous-volume vectorisé extrait (1D, de taille xr.size() * yr.size() * zr.size()).</returns>
std::vector<double> ExtractSubVolume(const std::vector<double>& data, const Vec3I& dim, const std::vector<size_t>& xrange,
                                     const std::vector<size_t>& yrange, const std::vector<size_t>& zrange);

//======================================================
//==================== Manipulation ====================
//======================================================

/// <summary>Calcule le centre de masse (x, y) d'une image 2D carrée de taille (nxn).</summary>
/// <param name="data">Données des pixels (valeurs d’intensité).</param>
/// <param name="n">Largeur/hauteur de l'image carrée (nxn).</param>
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
/// <param name="n"> Taille de la matrice (nxn).</param>
void InverseGaussJordan(std::vector<double>& in, std::vector<double>& out, size_t n);

/// <summary>Inversion de matrice (méthode LU (sans pivot) + systèmes triangulaires).</summary>
/// <param name="in">Matrice à inverser (format vectorisé ligne par ligne, modifiée en place).</param>
/// <param name="out">Matrice résultat de l'inversion (même format).</param>
/// <param name="n"> Taille de la matrice (nxn).</param>
void InverseLU(std::vector<double>& in, std::vector<double>& out, size_t n);

/// <summary>
/// Génère un vecteur d’indices réfléchis de longueur <c>k0</c>, utilisé pour initialiser les filtres exponentiels symétriques.
///
/// Le motif généré est basé sur la réflexion de la séquence [0, 1, ..., n-1] : [0, 1, ..., n-1, n-2, ..., 1], répété autant de fois que nécessaire.
///
/// <para>Exemple : ReflectionIndices(5, 10) ⇒ [0,1,2,3,4,3,2,1,0,1]</para>
///
/// </summary>
/// <param name="n">Taille de la dimension originale à réfléchir.</param>
/// <param name="k0">Longueur finale désirée du vecteur réfléchi.</param>
/// <returns>Vecteur d’indices réfléchis de taille <c>k0</c>.</returns>
std::vector<size_t> ReflectionIndices(size_t n, size_t k0);

/// <summary>
/// Génère un vecteur d’indices réfléchis (type miroir) pour le padding.
/// Par exemple pour N = 5 et pad = 1 : retourne [1,0,1,2,3,4,3,2]
/// </summary>
/// <param name="n">Taille de la dimension originale.</param>
/// <param name="pad">Nombre d’éléments à ajouter de chaque côté.</param>
/// <returns>Vecteur d’indices réfléchis de taille N + 2*pad.</returns>
std::vector<size_t> ReflectIndex1D(size_t n, size_t pad);
} // namespace palm_tracer::matrix
