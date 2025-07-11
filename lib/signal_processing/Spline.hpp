///-------------------------------------------------------------------------------------------------
///
/// \file Spline.hpp
/// \brief Classe de calibration de la spline pour l'ajustement.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
/// \remarks Ce code est basé sur celui de SMAP <a href="https://github.com/jries/SMAP">github.com/jries/SMAP</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <string>
#include <vector>

#include "Structs.hpp"

namespace palm_tracer {
/// <summary>Classe de calibration de la spline pour l'ajustement.</summary>
class Spline {
public:
	size_t imSize_{},
	       stackSize_{},
	       beadImSize_{},
	       beadStackSize_{},
	       smallBeadImSize_{},
	       smallBeadStackSize_{},

	       distanceHalf_{},
	       dz_{},
	       smoothZ_{};

	Coord3D size_{}, beadSize_{}, midBead_{}, smallBeadSize_{};
	std::vector<size_t> smallXIndices_, smallZIndices_;

	Spline()  = default;
	~Spline() = default;

	void Calibrate(const uint16_t* stack, const std::vector<Coord>& beads, const Coord3D& size, size_t roi, size_t distance, size_t dz, size_t smoothZ);

	void Load(const std::string& filename);

	void Save(const std::string& filename);

	// Future protected
	void SetMembers(const Coord3D& size, size_t roi, size_t distance, size_t dz, size_t smoothZ);

	std::vector<std::vector<double>> ExtractBeads(const uint16_t* stack, const std::vector<Coord>& beads) const;

	[[nodiscard]] std::vector<double> GetSubVolume(const std::vector<double>& volume, bool smallZ) const;

	void RegisterPSF(std::vector<std::vector<double>>& beads);

	std::vector<Vec3> AlignToRef(const std::vector<double>& ref, std::vector<std::vector<double>>& beads, bool smallZ = false) const;

	std::vector<std::vector<double>> softNormalization(std::vector<std::vector<double>>& beads);

	void GetOverlap(const std::vector<std::vector<double>>& beads, const std::vector<Vec3>& shift, std::vector<bool>& good, std::vector<double>& weight,
	                double& norm, bool smallZ = false) const;

protected:
};
} // namespace palm_tracer
