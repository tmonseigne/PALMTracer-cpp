///-------------------------------------------------------------------------------------------------
///
/// \file Spline.hpp
/// <summary>Classe de calibration de spline pour PSF 3D à partir d'un stack d'acquisition.
/// </summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 02/07/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// <remarks>Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.</remarks>
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "BSpline3D.hpp"
#include "Structs.hpp"

namespace palm_tracer {
/// <summary>Classe de calibration de spline pour PSF 3D pour l'ajustement à partir d'un stack de billes.</summary>
class Spline {
public:
	Vec3I dim_{},          ///< Dimensions complètes du volume 3D d'origine (X, Y, Z).
	      beadDim_{},      ///< Dimensions d'une ROI centrée sur une bille (cube impair).
	      midBead_{},      ///< Coordonnées du centre de la ROI (en voxels).
	      smallBeadDim_{}; ///< Dimensions de la zone centrale de la ROI utilisée pour l'alignement/normalisation.

	size_t imSize_{},          ///< Taille d'une image 2D (X*Y).
	       stackSize_{},       ///< Taille totale du stack (X*Y*Z).
	       beadImSize_{},      ///< Taille d'une coupe de la ROI (X*Y).
	       beadSize_{},        ///< Taille totale d'une ROI (X*Y*Z).
	       smallBeadImSize_{}, ///< Taille d'une coupe centrale de ROI.
	       smallBeadSize_{},   ///< Taille totale de la zone centrale de ROI.

	       distanceHalf_{}, ///< Moitié de la fenêtre autour du centre pour l'alignement Z.
	       dz_{};           ///< Pas d'acquisition sur Z (en nm).

	std::vector<size_t> smallXIndices_, ///< Indices X de la zone centrale utilisée pour la comparaison.
	                    smallZIndices_, ///< Indices Z de la zone centrale utilisée pour l'alignement.
	                    zIndices_;      ///< Indices Z complets de la ROI.

	Spline()  = default; ///< Constructeur de la classe <see cref="Spline"/>.
	~Spline() = default; ///< Destructeur de la classe <see cref="Spline"/>.

	/// <summary>Calibre une PSF 3D spline à partir d'une pile d'images et positions de billes.</summary>
	/// <param name="stack">Volume brut d'acquisition (stack 3D d'uint16_t).</param>
	/// <param name="beads">Liste des positions XY des billes.</param>
	/// <param name="dim">Dimensions complètes de la pile (X, Y, Z).</param>
	/// <param name="roi">Taille carrée de la ROI autour de chaque bille (sera rendue impaire automatiquement).</param>
	/// <param name="distance">Nombre de plans Z utilisés autour du centre (pour smallZ).</param>
	/// <param name="dz">Pas Z en nm (utilisé pour pondérer la régularisation de la spline sur Z).</param>
	/// <remarks>Cette fonction calcule une PSF moyenne spline à partir d'un ensemble de billes, aligne, filtre et normalise.</remarks>
	void Calibrate(const uint16_t* stack, const std::vector<Vec2I>& beads, const Vec3I& dim, size_t roi, size_t distance, size_t dz);

	/// <summary>Charge une spline depuis un fichier (non implémenté ici).</summary>
	void Load(const std::string& filename) = delete;

	/// <summary>Sauvegarde la spline dans un fichier (non implémenté ici).</summary>
	void Save(const std::string& filename) = delete;

	/// <summary>Retourne les coefficients finaux de la spline.</summary>
	/// <returns>Vecteur 1D des coefficients (dimFinal x 64).</returns>
	[[nodiscard]] std::vector<double> GetCoeffs() const { return bs_.coeffsFinal_; }

	// Future protected

	/// <summary>Initialise les paramètres internes (dimensions, indices, etc).</summary>
	/// <param name="dim">Dimensions complètes de la pile (X, Y, Z).</param>
	/// <param name="roi">Taille carrée de la ROI autour de chaque bille (sera rendue impaire automatiquement).</param>
	/// <param name="distance">Nombre de plans Z utilisés autour du centre (pour smallZ).</param>
	/// <param name="dz">Pas Z en nm (utilisé pour pondérer la régularisation de la spline sur Z).</param>
	void SetMembers(const Vec3I& dim, size_t roi, size_t distance, size_t dz);

	/// <summary>Extrait les volumes centrés autour de chaque bille depuis le stack brut.</summary>
	/// <param name="stack">Volume brut d'acquisition (stack 3D d'uint16_t).</param>
	/// <param name="beads">Liste des positions XY des billes.</param>
	/// <returns>Vecteur de volumes double (un par bille, linéarisé).</returns>
	std::vector<std::vector<double>> ExtractBeads(const uint16_t* stack, const std::vector<Vec2I>& beads) const;

	/// <summary>Calcule la PSF moyenne alignée à partir d'une liste de sous-volumes.</summary>
	/// <param name="beads">Volumes individuels des billes.</param>
	/// <returns>Volume moyen des billes.</returns>
	std::vector<double> GetPSF(std::vector<std::vector<double>>& beads) const;

	/// <summary>Aligne les billes sur la PSF de référence par shift subpixel.</summary>
	/// <param name="ref">Volume de référence (moyenne des billes).</param>
	/// <param name="beads">Volumes des billes à aligner (modifiés en place).</param>
	/// <param name="smallZ">Si vrai, utilise un sous-volume centré en Z.</param>
	/// <returns>Liste des vecteurs de décalage appliqués.</returns>
	std::vector<Vec3D> Align(const std::vector<double>& ref, std::vector<std::vector<double>>& beads, bool smallZ = false) const;

	/// <summary>Normalisation douce des intensités des billes, par rapport à la PSF moyenne locale.</summary>
	/// <param name="beads">Volumes individuels.</param>
	/// <returns>Volumes normalisés.</returns>
	[[nodiscard]] std::vector<std::vector<double>> SoftNormalization(const std::vector<std::vector<double>>& beads) const;

	/// <summary>Met à jour les poids des billes selon leur RMS d'erreur vis-à-vis de la référence.</summary>
	/// <param name="beads">Volumes (normalisés) des billes.</param>
	/// <param name="shift">Décalages appliqués pour chaque bille.</param>
	/// <param name="weights">Poids mis à jour (NaN si rejeté).</param>
	/// <param name="norm">Valeur de normalisation globale.</param>
	/// <param name="smallZ">Si vrai : réduction sur Z avant comparaison.</param>
	void GetOverlap(const std::vector<std::vector<double>>& beads, const std::vector<Vec3D>& shift, std::vector<double>& weights,
	                double& norm, bool smallZ = false) const;

protected:
	BSpline3D bs_; ///< Spline 3D finale obtenue après calibration.
};
} // namespace palm_tracer
