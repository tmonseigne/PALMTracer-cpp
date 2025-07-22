///-------------------------------------------------------------------------------------------------
///
/// \file Alignment.hpp
/// \brief Classe pour l'alignement optique d'une pile d'images 3D avec déformation non linéaire.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
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
	/// <summary>Constructeur de la classe <c>Alignment</c>.</summary>
	Alignment() = default;
	/// <summary>Destructeur de la classe <c>Alignment</c>.</summary>
	~Alignment() = default;

	/// <summary>Mode d'alignement avec calcul des coefficients à partir d'images de référence.</summary>
	/// <param name="highSNR">Image de référence (alignée) avec un haut rapport signal/bruit.</param>
	/// <param name="lowSNR">Image de référence (désalignée) avec un haut faible signal/bruit.</param>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la zone mémoire de sortie (même taille que `input`).</param>
	/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
	/// <param name="upsampling">Taux de suréchantillonnage à appliquer avant recalage.</param>
	/// <remarks>Cette méthode lance une exception, car le calcul automatique n'est pas implémenté dans cette version.</remarks>
	void Align(const uint16_t* highSNR, const uint16_t* lowSNR, const uint16_t* input, uint16_t* output, const Coord3D& dim, size_t upsampling = 1);

	/// <summary>Mode d'alignement avec coefficients déjà calculés.</summary>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la zone mémoire de sortie (même taille que `input`).</param>
	/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
	/// <param name="factors">Coefficients du modèle non linéaire, séparés pour X et Y (2x10 coefficients).</param>
	/// <param name="upsampling">Taux de suréchantillonnage à appliquer avant recalage (>=1).</param>
	void Align(const uint16_t* input, uint16_t* output, const Coord3D& dim, const std::array<std::array<double, 10>, 2>& factors, size_t upsampling = 1);

protected:
	/// <summary>Structure décrivant la contribution d'un pixel source à un pixel cible.</summary>
	struct Contribution {
		size_t index;  ///< Index dans l'image source (flattened).
		double weight; ///< Poids associé au pixel source.
	};

	const uint16_t *highSNR_{},                           ///< Image de référence alignée (fort SNR).
	               *lowSNR_{},                            ///< Image de référence désalignée (faible SNR).
	               *input_{};                             ///< Pointeur vers la pile d'entrée à aligner.
	uint16_t* output_{};                                  ///< Pointeur vers la zone de sortie (résultat de l'alignement).
	Coord3D dim_{};                                       ///< Dimensions du volume 3D (width=X, height=Y, depth=Z).
	size_t sizeIm_{},                                     ///< Nombre total de pixels d'une image (height × width).
	       sizeStack_{},                                  ///< Nombre total de pixels (height × width × planes).
	       upsampling_{},                                 ///< Taux de suréchantillonnage.
	       upsampling2_{};                                ///< Taux de suréchantillonnage au carré.
	std::array<std::array<double, 10>, 2> factors_{};     ///< Coefficients du modèle de déformation polynomiale.
	std::vector<std::vector<Contribution>> transformMap_; ///< Carte de transformation pré-calculée (chaque pixel cible contiendra une liste de contributions).

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
	/// Cette carte contient, pour chaque pixel de l’image cible, les coordonnées sources correspondantes ainsi que les poids associés à l’interpolation.</summary>
	/// <remarks>Utilise un splatting inverse avec suréchantillonnage.</remarks>
	void precomputeTransformMap();

	/// <summary>Applique les contributions pré-calculées sur toutes les images du volume.</summary>
	/// <remarks>Utilise OpenMP pour paralléliser le traitement sur chaque pixel de sortie.</remarks>
	void applyTransform() const;

	/// <summary>Réalise une découpe polygonale (algorithme de Sutherland-Hodgman) pour estimer l'aire de recouvrement.</summary>
	/// <param name="poly">Polygone d'origine.</param>
	/// <param name="clipped">Polygone résultant dans le pixel cible.</param>
	/// <param name="xmin">Coordonnée X du coin supérieur gauche du pixel cible.</param>
	/// <param name="ymin">Coordonnée Y du coin supérieur gauche du pixel cible.</param>
	/// <param name="pixelSize">Taille d'un pixel (souvent 1 / suréchantillonnage).</param>
	static void clipPolygon(const std::vector<Vec2>& poly, std::vector<Vec2>& clipped, double xmin, double ymin, double pixelSize);

	/// <summary>Calcule l'aire d'un polygone (utile pour les poids de contribution).</summary>
	/// <param name="poly">Polygone (doit être convexe pour validité des poids).</param>
	/// <returns>Aire (toujours positive, avec normalisation par le suréchantillonnage).</returns>
	[[nodiscard]] double getArea(const std::vector<Vec2>& poly) const;
};
} // namespace palm_tracer
