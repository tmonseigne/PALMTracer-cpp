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
	size_t roi_{},
	       roiHalf_{},
	       roi2_{},
	       height_{},
	       width_{},
	       planes_{},
	       planesHalf_{},
	       sizeIm_{},
	       sizeStack_{},
	       sizeBeads_{},
	       distance_{},
	       distanceHalf_{},
	       dz_{},
	       smoothZ_{},
	       zStart_,
	       zEnd_,
	       zN_;
	std::vector<size_t> roiCenterIndices;

	Spline()  = default;
	~Spline() = default;

	void Calibrate(const uint16_t* stack, const std::vector<Coord>& beads, size_t height, size_t width, size_t planes, size_t roi,
	               size_t distance, size_t dz, size_t smoothZ);

	void Load(const std::string& filename);

	void Save(const std::string& filename);

protected:
	void setMembers(size_t height, size_t width, size_t planes, size_t roi, size_t distance, size_t dz, size_t smoothZ);

	std::vector<std::vector<double>> extractBeads(const uint16_t* stack, const std::vector<Coord>& beads) const;

	void registerPSF(std::vector<std::vector<double>>& beads) const;

	void get3DShift(const std::vector<double>& ref, const std::vector<double>& target, Point3D& shift) const;

	static void idToPoint3D(size_t id, Point3D& pos, size_t x, size_t y);
};
} // namespace palm_tracer
