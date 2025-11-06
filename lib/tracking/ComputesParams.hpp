///-------------------------------------------------------------------------------------------------
///
/// \file ComputesParams.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <iostream>

namespace palm_tracer {
/// <summary>Modes de croissance de la distance maximale de recherche lors de la reconnexion de trajectoires séparées par scintillement.</summary>
enum class blinking_mode : uint8_t {
	immobile,  ///< Distance maximale constante avec le temps (recherche dans un rayon fixe).
	diffusive, ///< Distance maximale croissante en racine carrée du temps écoulé.
	linear     ///< Distance maximale croissante linéairement avec le temps écoulé.
};

/// <summary>Mode du modèle d'ajustement pour les calculs statistiques des trajectoires.
/// L'ajustement cherche le vecteur de paramètres <c>p</c> qui minimise :
/// \f[ \Phi(p) = \sum_{i=0}^{N-1} \big(y_i - f(x_i; p)\big)^2 \f]
/// (moindres carrés non pondérés). La qualité peut être évaluée via <c>computeR2()</c>.
/// </summary>
/// <remarks>Plus de cas avec le domaine spectral et les histogrammes seront implémantés plus tard :
/// Spectral : Bi-dimensional Polynomial, Translation, Affine, Affine (Translation), Polywarp.
/// Histogramm : Gaussian, Bi-Gaussian (donc Histo_Gaussian, Histo_BiGaussian).
///	Autres méthodes appelées directement : Third Order Polynomial, 3D Astygmatism, NeNa, Gaussian.
/// </remarks>
enum class fit_mode : uint8_t {
	linear,     ///< \f$p = [a, b]\f$ pour optimiser la fonction \f$ f(x) = a \times x + b\f$
	power,      ///< \f$p = [\alpha, a, b]\f$ pour optimiser la fonction \f$ f(x) = a \times x^{α} + b\f$
	exponential ///< \f$p = [x_0, a, \tau, b]\f$ pour optimiser la fonction \f$ f(x) = (b - a) \times \exp\Big(-\tfrac{x - x_0}{\tau}\Big) + a \f$
};

/// <summary>Structure contenant les paramètres de calcul pour le module de tracking.</summary>
struct ComputesParams {
	fit_mode mode       = fit_mode::linear; ///< Mode d'ajustement utilisé pour les calculs statistiques (voir <see cref="fit_mode"/>).
	size_t fitLength    = 0;                ///< Longueur de la fenêtre de calcul (nombre d'échantillons).
	bool is3D           = false;            ///< Effectuer des calculs en 3D (true) ou en 2D (false).
	double pixelSize    = 0.160;            ///< Taille d'un pixel en unitées réelles (micromètre).
	double exposureTime = 0.05;             ///< Temps d'exposition en secondes.
};

/// <summary>Surcharge de l'opérateur d'affichage pour la structure ComputesParams. Format : (...) </summary>
inline std::ostream& operator<<(std::ostream& os, const ComputesParams& p)
{
	os << "(" << static_cast<int>(p.mode) << ", " << p.fitLength << ", " << p.is3D << ", " << p.pixelSize << ", " << p.exposureTime << ")";
	return os;
}
} // namespace palm_tracer
