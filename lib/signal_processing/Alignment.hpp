///-------------------------------------------------------------------------------------------------
///
/// \file Alignment.hpp
/// <summary>Classe pour l'alignement optique d'une pile d'images 3D avec déformation non linéaire.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 02/06/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <array>
#include <cstdint>
#include <vector>
#include "Structs.hpp"

namespace palm_tracer {
/// <summary>Classe pour l'alignement optique d'une pile d'images 3D avec déformation non linéaire.</summary>
///
/// Cette classe applique un recalage non linéaire à une pile d'images 3D (Z-stacks) en utilisant un modèle polynomial 3e degré.
/// Deux modes d'utilisation sont possibles :
/// - soit le calcul des coefficients de transformation à partir d'une image référence haute et basse SNR,
/// - soit l'application directe de coefficients pré-calculés.
class Alignment {
public:
	/// <summary>Constructeur de la classe <see cref="Alignment"/>.</summary>
	Alignment() = default;
	/// <summary>Destructeur de la classe <see cref="Alignment"/>.</summary>
	~Alignment() = default;

	/// <summary>Mode d'alignement avec calcul des coefficients à partir d'images de référence.</summary>
	/// <param name="highSNR">Image de référence (alignée) avec un haut rapport signal/bruit.</param>
	/// <param name="lowSNR">Image de référence (désalignée) avec un faible rapport signal/bruit.</param>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la sortie, de taille <c>dim.GetSize() * upsampling * upsampling</c>.</param>
	/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
	/// <param name="upsampling">Taux de suréchantillonnage à appliquer avant recalage.</param>
	/// <remarks>Cette méthode lance une exception, car le calcul automatique n'est pas implémenté dans cette version.</remarks>
	void Align(const uint16_t* highSNR, const uint16_t* lowSNR, const uint16_t* input, uint16_t* output, const Vec3I& dim, size_t upsampling = 1);

	/// <summary>Mode d'alignement avec coefficients déjà calculés.</summary>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la sortie, de taille <c>dim.GetSize() * upsampling * upsampling</c>.</param>
	/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
	/// <param name="factors">Coefficients du modèle non linéaire, séparés pour X et Y (2x10 coefficients).</param>
	/// <param name="upsampling">Taux de suréchantillonnage à appliquer avant recalage (>=1).</param>
	void Align(const uint16_t* input, uint16_t* output, const Vec3I& dim, const std::array<std::array<double, 10>, 2>& factors, size_t upsampling = 1);

protected:
	/// <summary>Structure décrivant la contribution d'un pixel source à un pixel cible.</summary>
	struct Contribution {
		size_t index;  ///< Index dans l'image source (flattened).
		double weight; ///< Poids associé au pixel source.
	};

	const uint16_t *highSNR_{},                       ///< Image de référence alignée (fort SNR).
	               *lowSNR_{},                        ///< Image de référence désalignée (faible SNR).
	               *input_{};                         ///< Pointeur vers la pile d'entrée à aligner.
	uint16_t* output_{};                              ///< Pointeur vers la zone de sortie (résultat de l'alignement).
	Vec3I dim_{};                                     ///< Dimensions du volume 3D (width=X, height=Y, depth=Z).
	size_t sizeIm_{},                                 ///< Nombre total de pixels d'une image (height × width).
	       sizeStack_{},                              ///< Nombre total de pixels (height × width × planes).
	       upsampling_{},                             ///< Taux de suréchantillonnage.
	       upsampling2_{};                            ///< Taux de suréchantillonnage au carré.
	std::array<std::array<double, 10>, 2> factors_{}; ///< Coefficients du modèle de déformation polynomiale.
	std::vector<size_t> offsets_;                     ///< Débuts des listes CSR, suivis du nombre total de contributions.
	std::vector<Contribution> contributions_;         ///< Contributions contiguës, triées par indice source pour chaque cible.

	/// <summary>Initialise la mémoire de sortie à 0.</summary>
	void initOutput() const { for (size_t i = 0; i < sizeStack_; ++i) { output_[i] = 0; } }

	/// <summary>Méthode statique pour calculer les coefficients de recalage (non implémentée ici).</summary>
	/// <remarks>Lance une exception. Cette méthode est un placeholder pour une version future.</remarks>
	static void computeFactor();

	/// <summary>Applique le modèle de déformation non linéaire aux coordonnées (x, y).</summary>
	/// <param name="oldX">Coordonnée X d'origine.</param>
	/// <param name="oldY">Coordonnée Y d'origine.</param>
	/// <param name="newX">Coordonnée X recalculée.</param>
	/// <param name="newY">Coordonnée Y recalculée.</param>
	void computeCoord(double oldX, double oldY, double& newX, double& newY) const;

	/// <summary>Applique l'alignement sur la pile d'images 3D.</summary>
	/// <remarks>Cette méthode prépare la carte de transformation et applique la réallocation.</remarks>
	void align();

	/// <summary>>Pré-calcul des contributions pour chaque pixel cible en fonction de la déformation pour une interpolation plus rapide.
	///
	/// Cette carte contient, pour chaque pixel de l'image cible, les coordonnées sources correspondantes ainsi que les poids associés à
	/// l'interpolation.</summary>
	/// <remarks>Les contributions sont sommées en double par indice source croissant, indépendamment des threads.
	/// Les coefficients non finis donnent une carte vide. Les cellules à coins non finis, dégénérées ou non convexes sont ignorées.
	/// Les réflexions sont acceptées. Aucun inverse global n'est nécessaire : les recouvrements sont additionnés.
	/// La carte CSR est construite par lignes sources indépendantes, puis comptage, préfixe et remplissage stable.</remarks>
	void precomputeTransformMap();

	/// <summary>Applique les contributions pré-calculées sur toutes les images du volume.</summary>
	/// <remarks>Arrondi au plus proche, demi-entiers vers le haut, puis saturation dans <c>[0, 65535]</c>.
	/// Les sommes non positives ou NaN donnent zéro ; les débordements positifs donnent 65535.
	/// Utilise OpenMP pour paralléliser le traitement sur chaque pixel de sortie.</remarks>
	void applyTransform() const;

	/// <summary>Réalise une découpe polygonale (algorithme de Sutherland-Hodgman) pour estimer l'aire de recouvrement.</summary>
	/// <param name="poly">Polygone d'origine.</param>
	/// <param name="clipped">Polygone résultant dans le pixel cible.</param>
	/// <param name="xmin">Coordonnée X du coin supérieur gauche du pixel cible.</param>
	/// <param name="ymin">Coordonnée Y du coin supérieur gauche du pixel cible.</param>
	/// <param name="pixelSize">Taille d'un pixel (souvent 1 / suréchantillonnage).</param>
	static void clipPolygon(const std::vector<Vec2D>& poly, std::vector<Vec2D>& clipped, double xmin, double ymin, double pixelSize);

	/// <summary>Calcule l'aire d'un polygone (utile pour les poids de contribution).</summary>
	/// <param name="poly">Polygone (doit être convexe pour validité des poids).</param>
	/// <returns>Aire (toujours positive, avec normalisation par le suréchantillonnage).</returns>
	[[nodiscard]] double getArea(const std::vector<Vec2D>& poly) const;
};
} // namespace palm_tracer
