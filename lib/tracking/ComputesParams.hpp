///-------------------------------------------------------------------------------------------------
///
/// \file ComputesParams.hpp
/// <summary>Définition des paramètres de calcul utilisés par le module de suivi et d'analyse des trajectoires.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <cstdint>
#include <iostream>

namespace palm_tracer::tracking {
/// <summary>Modes de croissance de la distance maximale de recherche lors de la reconnexion de trajectoires séparées par scintillement.</summary>
enum class blinking_mode : uint8_t {
	immobile,  ///< Distance maximale constante avec le temps (recherche dans un rayon fixe).
	diffusive, ///< Distance maximale croissant en racine carrée du temps écoulé.
	linear     ///< Distance maximale croissant linéairement avec le temps écoulé.
};

/// <summary>Mode du modèle d'ajustement pour l'analyse statistique des trajectoires.
/// L'ajustement cherche le vecteur de paramètres <c>p</c> qui minimise la somme des carrés :
/// \f[ \Phi(p) = \sum_{i=0}^{N-1}\big(y_i - f(x_i; p)\big)^2 \f] (moindres carrés non pondérés).
/// La qualité peut être évaluée via <see cref="Fit::computeR2"/>.</summary>
/// <remarks>
/// D'autres cas (espace spectral, histogrammes) pourront être ajoutés ultérieurement :
/// - spectral : polynôme bidimensionnel, translation, affine (avec/sans translation), polywarp ;
/// - histogrammes : gaussien, bi-gaussien ; et méthodes directes diverses (polynôme d'ordre 3, astigmatisme 3D, NeNa, gaussien).
/// </remarks>
enum class fit_mode : uint8_t {
	linear,     ///< Paramètres \f$ p = [a,\, b] \f$ pour \f$ f(x) = a\,x + b \f$.
	power,      ///< Paramètres \f$ p = [\alpha,\, a,\, b] \f$ pour \f$ f(x) = a\,x^{\alpha} + b \f$.
	exponential ///< Paramètres \f$ p = [x_0,\, a,\, \tau,\, b] \f$ pour \f$ f(x) = (b - a)\,\exp\!\big(-(x - x_0)/\tau\big) + a \f$.
};

/// <summary>Structure contenant les paramètres de calcul pour le module de tracking.</summary>
struct ComputesParams {
	fit_mode mode       = fit_mode::linear; ///< Mode d'ajustement utilisé pour les calculs statistiques (voir <see cref="fit_mode"/>).
	size_t fitLength    = 0;                ///< Longueur de la fenêtre de calcul (nombre d'échantillons).
	bool is3D           = false;            ///< Active le calcul 3D si <c>true</c> (sinon 2D).
	double pixelSize    = 0.160;            ///< Taille d'un pixel en unités réelles (micromètre).
	double exposureTime = 0.05;             ///< Temps d'exposition (seconde), noté \f$ \Delta t \f$.
};

/// <summary>Surcharge de l'opérateur d'affichage pour <see cref="ComputesParams"/>. Format : <c>(mode, fitLength, is3D, pixelSize, exposureTime)</c>.</summary>
inline std::ostream& operator<<(std::ostream& os, const ComputesParams& p)
{
	os << "(" << static_cast<int>(p.mode) << ", " << p.fitLength << ", " << p.is3D << ", " << p.pixelSize << ", " << p.exposureTime << ")";
	return os;
}
} // namespace palm_tracer::tracking
