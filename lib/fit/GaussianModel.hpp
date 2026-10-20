#pragma once

#include <span>

namespace palm_tracer::fit {
/// <summary>État gaussien immuable préparé pour plusieurs pixels d'un même jeu de paramètres.</summary>
/// <remarks>Les largeurs anisotropes suivent la convention interne sans facteur 2 dans l'exponentielle.
/// Ce modèle est préparé une fois par état de l’optimiseur puis évalué sur tous les pixels de la ROI.</remarks>
class GaussianModel {
public:
	/// <summary>Prépare les constantes d'une gaussienne isotrope.</summary>
	/// <param name="offset">Fond du modèle.</param>
	/// <param name="amplitude">Amplitude crête.</param>
	/// <param name="sigma">Écart-type en pixels.</param>
	/// <param name="hasSigma">Active la dérivée de sigma en cinquième position.</param>
	GaussianModel(double offset, double amplitude, double sigma, bool hasSigma);

	/// <summary>Prépare les constantes d'une gaussienne anisotrope orientée.</summary>
	/// <param name="offset">Fond du modèle.</param>
	/// <param name="amplitude">Amplitude crête.</param>
	/// <param name="sigmaX">Largeur interne horizontale.</param>
	/// <param name="sigmaY">Largeur interne verticale.</param>
	/// <param name="theta">Angle en radians.</param>
	/// <param name="hasTheta">Active la dérivée de l'angle en septième position.</param>
	GaussianModel(double offset, double amplitude, double sigmaX, double sigmaY, double theta, bool hasTheta);

	/// <summary>Évalue la valeur et, facultativement, son gradient en partageant la même exponentielle.</summary>
	/// <param name="x">Coordonnée horizontale relative au centre.</param>
	/// <param name="y">Coordonnée verticale relative au centre.</param>
	/// <param name="gradient">Vue vide pour la valeur seule, sinon au moins 4/5 éléments isotropes ou 6/7 anisotropes.
	/// Les éléments supplémentaires restent inchangés.</param>
	/// <returns>Valeur de la PSF au pixel.</returns>
	/// <exception cref="std::invalid_argument">La vue non vide ne peut pas contenir le gradient demandé.</exception>
	[[nodiscard]] double Evaluate(double x, double y, std::span<double> gradient = {}) const;

private:
	const double offset_, amplitude_;                                        ///< Fond et amplitude de l'état courant.
	const double sigmaX_, sigmaY_, sigmaX2_, sigmaY2_;                       ///< Largeurs et carrés préparés.
	const double cosTheta_, sinTheta_, cosTheta2_, sinTheta2_, cosSinTheta_; ///< Constantes de rotation.
	const bool isotropic_, hasSigma_, hasTheta_;                             ///< Modèle et composantes libres du gradient.
};
} // namespace palm_tracer::fit
