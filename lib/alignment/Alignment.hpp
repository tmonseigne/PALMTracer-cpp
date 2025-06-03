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
	/// <param name="upscale">Taux de suréchantillonnage à appliquer avant recalage.</param>
	/// <param name="upscaleMode">Méthode de suréchantillonnage (utilisable avec OpenCV ou autre).</param>
	/// <param name="upscaleSource">Si vrai, les coefficients sont calculés sur les images de référence suréchantillonnées.</param>
	/// <param name="spread">Si vrai, active la dispersion sub-pixellaire (ex. bilinéaire).</param>
	void Align(const unsigned short* highSNR, const unsigned short* lowSNR, const unsigned short* input, unsigned short* output,
	           size_t height, size_t width, size_t planes, size_t upscale = 1, size_t upscaleMode = 0, bool upscaleSource = true, bool spread = true);

	/// <summary>Constructeur alternatif de la classe <c>Fitting</c>. Applique directement des coefficients d'alignement fournis.</summary>
	/// <param name="input">Données d'entrée (pile d'images à réaligner).</param>
	/// <param name="output">Pointeur vers la zone mémoire de sortie (même taille que `input`).</param>
	/// <param name="height">Hauteur de chaque plan (en pixels).</param>
	/// <param name="width">Largeur de chaque plan (en pixels).</param>
	/// <param name="planes">Nombre de plans (profondeur de la pile).</param>
	/// <param name="factors">Coefficients du modèle non linéaire (degrés 1 à 3).</param>
	/// <param name="upscale">Taux de suréchantillonnage à appliquer avant recalage.</param>
	/// <param name="upscaleMode">Méthode de suréchantillonnage (utilisable avec OpenCV ou autre).</param>
	/// <param name="upscaleSource">Si vrai, les coefficients sont calculés sur les images de référence suréchantillonnées.</param>
	/// <param name="spread">Si vrai, active la dispersion sub-pixellaire (ex. bilinéaire).</param>
	void Align(const unsigned short* input, unsigned short* output, size_t height, size_t width,
	           size_t planes, const std::array<std::array<double, 10>, 2>& factors,
	           size_t upscale = 1, size_t upscaleMode = 0, bool upscaleSource = true, bool spread = true);

protected:
	const unsigned short *highSNR_,                 ///< Image de référence alignée (fort SNR)
	                     *lowSNR_,                  ///< Image de référence désalignée (faible SNR)
	                     *input_;                   ///< Données d'entrée (pile à aligner)
	unsigned short* output_;                        ///< Résultat de l'alignement
	size_t height_,                                 ///< Hauteur des images (pixels)
	       width_,                                  ///< Largeur des images (pixels)
	       planes_,                                 ///< Nombre de plans dans la pile
	       size_,                                   ///< Nombre total de pixels (height × width × planes)
	       upscale_,                                ///< Taux de suréchantillonnage
	       upscaleMode_,                            ///< Méthode de suréchantillonnage
	       upsacleSource_,                          ///< Indique si les références ont été upscalées
	       spread_;                                 ///< Active la dispersion lors du rééchantillonnage
	std::array<std::array<double, 10>, 2> factors_; ///< Coefficients du modèle polynomial

	/// <summary>Initialise la mémoire de sortie à 0.</summary>
	void initOutput() const;

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

	/// <summary>Convertit des coordonnées flottantes (x, y) en index linéaire dans le plan.</summary>
	/// <param name="x">Coordonnée X (flottante).</param>
	/// <param name="y">Coordonnée Y (flottante).</param>
	/// <returns>Index linéaire dans l'image 2D.</returns>
	[[nodiscard]] size_t getCoord(double x, double y) const;

	/// <summary>Applique l'alignement sur la pile d'images 3D.</summary>
	///
	/// Redirige chaque pixel selon la transformation non linéaire.
	void align() const;
};
}
