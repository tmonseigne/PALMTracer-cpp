///-------------------------------------------------------------------------------------------------
/// \file FrameProcessor.hpp
/// <summary>Définition du pipeline réutilisable de traitement d’une image.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "ATrous.hpp"
#include "Segmentation2D.hpp"
#include "../Constants.hpp"
#include "../Structs.hpp"

namespace palm_tracer::localization {
/// <summary>Localisation interne, avec les colonnes historiques de sortie.</summary>
using LocalizationRow = std::array<double, localization_columns::Count>;
/// <summary>Localisations d'une image, triées par Y puis X après traitement.</summary>
using FrameLocalizations = std::vector<LocalizationRow>;

/// <summary>Enchaîne ondelettes, segmentation et ajustement pour une image.</summary>
/// <remarks>Réutilise ses espaces de travail pour des images de mêmes dimensions.
/// Une instance appartient exclusivement à un thread à la fois ; elle ne conserve pas les pixels des images traitées.</remarks>
class FrameProcessor {
public:
	/// <summary>Configure le traitement et alloue les espaces de travail.</summary>
	/// <param name="dim">Dimensions communes aux images.</param>
	/// <param name="maxPoints">Nombre maximal de localisations retenues par image.</param>
	/// <param name="watershedRatio">Seuil relatif de fusion des bassins.</param>
	/// <param name="fit">Mode historique : 1 à 4 pour les gaussiennes, 5 pour les splines, sinon sans ajustement.</param>
	/// <param name="fitParams">Paramètres empruntés, valides pendant la durée de vie du processeur.</param>
	FrameProcessor(const Vec2I& dim, size_t maxPoints, double watershedRatio,
	               size_t fit, const double* fitParams);

	/// <summary>Remplace les résultats par les localisations de l'image courante et les trie par Y puis X.</summary>
	/// <param name="image">Pixels contigus, au nombre de X × Y, valides pendant l'appel.</param>
	/// <param name="threshold">Seuil appliqué aux coefficients d'ondelettes.</param>
	/// <param name="points">Résultats limités à la capacité configurée, sans attribution des identifiants globaux.</param>
	void Process(const uint16_t* image, double threshold, FrameLocalizations& points);

private:
	Vec2I dim_;                   ///< Dimensions communes aux images.
	size_t maxPoints_;            ///< Capacité maximale par image.
	double watershedRatio_;       ///< Seuil relatif de fusion des bassins.
	size_t fit_;                  ///< Mode historique d'ajustement.
	const double* fitParams_;     ///< Paramètres empruntés pendant la durée de vie du processeur.
	ATrous decomposition_;        ///< Espace de travail des ondelettes.
	Segmentation2D segmentation_; ///< Espace de travail de segmentation.

	/// <summary>Ajuste les régions de la segmentation courante selon le mode configuré.</summary>
	/// <param name="image">Pixels de l'image ayant produit la segmentation.</param>
	/// <param name="points">Résultats remplacés par les localisations ajustées.</param>
	void fitRegions(const uint16_t* image, FrameLocalizations& points) const;
};
} // namespace palm_tracer::localization
