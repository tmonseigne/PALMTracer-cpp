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
#include <vector>

namespace palm_tracer {
/// <summary> Classe pour l'alignement optique d'une pile d'images 3D avec déformation non linéaire. </summary>
///
/// Cette classe applique un recalage non linéaire à une pile d'images 3D (Z-stacks) en utilisant un modèle polynomial 3e degré.
/// Elle peut soit estimer les coefficients à partir de deux images de référence (fort/faible SNR), soit appliquer des coefficients pré-calculés.
class Alignment {
public:
	/// <summary>Constructeur de la classe <c>Alignment</c>.</summary>
	Alignment() = default;
	/// <summary>Destructeur de la classe <c>Alignment</c>.</summary>
	~Alignment() = default;

	/// <summary>Constructeur de la classe <c>Fitting</c>. Calcule les coefficients d'alignement à partir de deux images de référence et lance l'alignement.</summary>
	/// <param name="highSNR">Image de référence (alignée) avec un haut rapport signal/bruit.</param>
	/// <param name="lowSNR">Image de référence (désalignée) avec un haut faible signal/bruit.</param>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la zone mémoire de sortie (même taille que `input`).</param>
	/// <param name="height">Hauteur de chaque plan (en pixels).</param>
	/// <param name="width">Largeur de chaque plan (en pixels).</param>
	/// <param name="planes">Nombre de plans (profondeur de la pile).</param>
	/// <param name="upsampling">Taux de suréchantillonnage à appliquer avant recalage.</param>
	void Align(const unsigned short* highSNR, const unsigned short* lowSNR, const unsigned short* input, unsigned short* output,
	           size_t height, size_t width, size_t planes, size_t upsampling = 1);

	/// <summary>Constructeur alternatif de la classe <c>Fitting</c>. Applique directement des coefficients d'alignement fournis.</summary>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la zone mémoire de sortie (même taille que `input`).</param>
	/// <param name="height">Hauteur de chaque plan (en pixels).</param>
	/// <param name="width">Largeur de chaque plan (en pixels).</param>
	/// <param name="planes">Nombre de plans (profondeur de la pile).</param>
	/// <param name="factors">Coefficients du modèle non linéaire (degrés 1 à 3).</param>
	/// <param name="upsampling">Taux de suréchantillonnage à appliquer avant recalage.</param>
	void Align(const unsigned short* input, unsigned short* output, size_t height, size_t width,
	           size_t planes, const std::array<std::array<double, 10>, 2>& factors, size_t upsampling = 1);

protected:
	/// <summary> Structure permettant de stocker les contributions sur la sortie. </summary>
	struct Contribution {
		size_t index;  ///< Index sur l'image en entrée.
		double weight; ///< Poids du pixel.
	};

	/// <summary> Structure permettant de stocker les coordonnées x et y en flottant.</summary>
	struct Coord {
		double x, y;
	};

	const unsigned short *highSNR_,                       ///< Image de référence alignée (fort SNR)
	                     *lowSNR_,                        ///< Image de référence désalignée (faible SNR)
	                     *input_;                         ///< Données d'entrée (pile à aligner)
	unsigned short* output_;                              ///< Résultat de l'alignement
	size_t height_,                                       ///< Hauteur des images (pixels)
	       width_,                                        ///< Largeur des images (pixels)
	       planes_,                                       ///< Nombre de plans dans la pile
	       sizeIm_,                                       ///< Nombre total de pixels d'une image (height × width)
	       size_,                                         ///< Nombre total de pixels (height × width × planes)
	       upsampling_;                                   ///< Taux de suréchantillonnage
	std::array<std::array<double, 10>, 2> factors_;       ///< Coefficients du modèle polynomial
	std::vector<std::vector<Contribution>> transformMap_; ///< chaque pixel cible contiendra une liste de contributions

	/// <summary>Initialise la mémoire de sortie à 0.</summary>
	void initOutput() const { for (size_t i = 0; i < size_; ++i) { output_[i] = 0; } }

	/// <summary>Fonction virtuelle de calcul des coefficients d'alignement. </summary>
	///
	/// Lance une exception, car non implémentée dans cette version.
	///
	static void computeFactor();

	/// <summary>Calcule les nouvelles coordonnées (x', y') selon le modèle non linéaire.</summary>
	/// <param name="oldX">Ancienne coordonnée X.</param>
	/// <param name="oldY">Ancienne coordonnée Y.</param>
	/// <param name="newX">Référence vers la nouvelle coordonnée X.</param>
	/// <param name="newY">Référence vers la nouvelle coordonnée Y.</param>
	void computeCoord(double oldX, double oldY, double& newX, double& newY) const;

	/// <summary>Applique l'alignement sur la pile d'images 3D.</summary>
	///
	/// Redirige chaque pixel selon la transformation non linéaire.
	void align();

	/// <summary>Pré-calcul d’une carte de transformation inverse pour une interpolation plus rapide.</summary>
	///
	/// Cette carte contient, pour chaque pixel de l’image cible, les coordonnées sources correspondantes
	/// ainsi que les poids associés à l’interpolation.
	void precomputeTransformMap();

	/// <summary>Applique la transformation pré-calculée à l’image d’entrée et produit l’image transformée.
	/// Cette fonction utilise la carte générée par `precomputeTransformMap`.</summary>
	void applyTransform() const;

	/// <summary>Coupe un polygone selon un carré de taille donnée.
	///
	/// Le carré est défini par son coin supérieur gauche `(xmin, ymin)` et sa taille `pixelSize`.
	/// Cette fonction utilise l’algorithme de Sutherland-Hodgman pour conserver la partie du polygone à l’intérieur du carré.</summary>
	/// <param name="poly">Polygone initial à découper.</param>
	/// <param name="clipped">Polygone résultant, après découpe.</param>
	/// <param name="xmin">Coordonnée X du coin supérieur gauche du carré.</param>
	/// <param name="ymin">Coordonnée Y du coin supérieur gauche du carré.</param>
	/// <param name="pixelSize">Taille du carré (largeur et hauteur identiques).</param>
	static void clipPolygon(const std::vector<Coord>& poly, std::vector<Coord>& clipped, double xmin, double ymin, double pixelSize);

	/// <summary>Calcule l’aire (signée) d’un polygone 2D.
	/// Le polygone peut être orienté : une aire négative indique un ordre antihoraire des sommets.</summary>
	/// <param name="poly">Liste des sommets du polygone (fermés ou non).</param>
	/// <returns>Aire du polygone en unités carrées.</returns>
	static double getArea(const std::vector<Coord>& poly);
};
}
