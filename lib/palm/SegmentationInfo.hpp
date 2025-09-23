///-------------------------------------------------------------------------------------------------
///
/// \file SegmentationInfo.hpp
/// \brief  Classe de lien pour traiter les données issues de la segmentation d'images 2D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

namespace palm_tracer {
/// <summary>Structure pour la fusion de régions (segmentation) avec potentiels et liens.</summary>
struct WSFusion {
	double maxPotentialWithNeighbor; ///< Potentiel maximal avec un voisin.
	double minPotential;             ///< Potentiel minimal.
	double maxPotential;             ///< Potentiel maximal.
	int x;                           ///< Position x.
	int y;                           ///< Position y.
	int structInContact;             ///< Index de la structure en contact.
	int structToBecome;              ///< Index de la structure à devenir.
	int newIndex;                    ///< Nouvel index assigné.
	bool locked;                     ///< Indique si verrouillé (pas modifiable).
};

/// <summary>Classe contenant des informations descriptives sur une région segmentée 2D.</summary>
class SegmentationInfo {
public:
	double x_,         ///< Position x.
	       y_,         ///< Position y.
	       z_,         ///< Position z.
	       intensity_, ///< Intensité totale dans la région.
	       surf_,      ///< Surface de la région.
	       mean_,      ///< Moyenne des intensités dans la région.
	       stdev_;     ///< Écart-type des intensités dans la région.

	/// <summary>Constructeur par défaut.</summary>
	SegmentationInfo() = default;

	/// <summary>Destructeur par défaut.</summary>
	~SegmentationInfo() = default;

	/// <summary>Constructeur de déplacement par défaut.</summary>
	SegmentationInfo(SegmentationInfo&&) noexcept = default;

	/// <summary>Opérateur d'affectation par déplacement par défaut.</summary>
	SegmentationInfo& operator=(SegmentationInfo&&) noexcept = default;

	/// <summary>Constructeur de copie par défaut.</summary>
	SegmentationInfo(const SegmentationInfo&) = default;

	/// <summary>Opérateur d'affectation par copie par défaut.</summary>
	SegmentationInfo& operator=(const SegmentationInfo&) = default;

	/// <summary>Réinitialise les données de la région.</summary>
	void Reset()
	{
		x_         = 0;
		y_         = 0;
		z_         = 0;
		surf_      = 0;
		stdev_     = 0;
		mean_      = 0;
		intensity_ = 0;
	}
};
} // namespace palm_tracer
