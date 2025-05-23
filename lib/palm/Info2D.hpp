///-------------------------------------------------------------------------------------------------
///
/// \file Info2D.hpp
/// \brief  Classe de lien pour traiter les données issues de la segmentation d'images 2D.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once

/// <summary>Structure représentant une coordonnée 2D en double précision.</summary>
struct Val2D {
	double x; ///< Coordonnée x.
	double y; ///< Coordonnée y.
};

/// <summary>Structure représentant une boîte englobante 2D par ses bornes.</summary>
struct Box2D {
	size_t xMin; ///< Borne minimale en x.
	size_t xMax; ///< Borne maximale en x.
	size_t yMin; ///< Borne minimale en y.
	size_t yMax; ///< Borne maximale en y.
};


/// <summary>Structure représentant un maximum local en 2D.</summary>
struct Max2D {
	double val; ///< Valeur du maximum.
	size_t x;   ///< Position x du maximum.
	size_t y;   ///< Position y du maximum.
};

/// <summary>Structure contenant les variances et covariance dans le plan 2D.</summary>
struct Variance2D {
	double xx; ///< Variance selon x.
	double yy; ///< Variance selon y.
	double xy; ///< Covariance entre x et y.
};

/// <summary>Structure pour la fusion de régions (segmentation) avec potentiels et liens.</summary>
struct WSFusion2D {
	double maxPotentialWithNeighbor; ///< Potentiel maximal avec un voisin.
	double minPotential;             ///< Potentiel minimal.
	double maxPotential;             ///< Potentiel maximal.
	int x;                           ///< Position x.
	int y;                           ///< Position y.
	int structInContact;             ///< Index de la structure en contact.
	int structToBecome;              ///< Index de la structure à devenir.
	int newIndex;                    ///< Nouvel index assigné.
	bool kept;                       ///< Indique si gardé après fusion.
	bool locked;                     ///< Indique si verrouillé (pas modifiable).
};

/// <summary>Classe contenant des informations descriptives sur une région segmentée 2D.</summary>
class Info2D {
public:
	size_t index_;      ///< Index identifiant la région.
	size_t peri_;       ///< Périmètre de la région.
	double surf_;       ///< Surface de la région.
	Max2D min_;         ///< Valeur minimale et sa position dans la région.
	Max2D max_;         ///< Valeur maximale et sa position dans la région.
	double stdev_;      ///< Écart-type des intensités dans la région.
	Val2D moment_;      ///< Moment géométrique (centre de masse).
	Val2D momentSurf_;  ///< Moment pondéré par la surface.
	double mean_;       ///< Moyenne des intensités dans la région.
	double intensity_;  ///< Intensité totale dans la région.
	Variance2D var_;    ///< Variances et covariance géométriques.
	Val2D sigma_;       ///< Écart-type géométrique.
	double angle_;      ///< Angle principal d'orientation.
	Box2D boundingBox_; ///< Boîte englobante de la région.

	/// <summary>Constructeur par défaut.</summary>
	Info2D() = default;

	/// <summary>Destructeur par défaut.</summary>
	~Info2D() = default;

	/// <summary>Constructeur de déplacement par défaut.</summary>
	Info2D(Info2D&&) noexcept = default;

	/// <summary>Opérateur d'affectation par déplacement par défaut.</summary>
	Info2D& operator=(Info2D&&) noexcept = default;

	/// <summary>Constructeur de copie par défaut.</summary>
	Info2D(const Info2D&) = default;

	/// <summary>Opérateur d'affectation par copie par défaut.</summary>
	Info2D& operator=(const Info2D&) = default;

	/// <summary> Réinitialise les données de la région.</summary>
	/// <param name="height">Hauteur maximale initiale (utilisée pour la boîte englobante)</param>
	/// <param name="width">Largeur maximale initiale (utilisée pour la boîte englobante)</param>
	void Reset(const size_t height = 0, const size_t width = 0)
	{
		peri_        = 0;
		surf_        = 0;
		min_         = {.val = 0, .x = 0, .y = 0};
		max_         = {.val = 0, .x = 0, .y = 0};
		stdev_       = 0;
		moment_      = {.x = 0, .y = 0};
		momentSurf_  = {.x = 0, .y = 0};
		mean_        = 0;
		intensity_   = 0;
		var_         = {.xx = 0, .yy = 0, .xy = 0};
		sigma_       = {.x = 0, .y = 0};
		angle_       = 0;
		boundingBox_ = {.xMin = height, .xMax = 0, .yMin = width, .yMax = 0};
	}
};
