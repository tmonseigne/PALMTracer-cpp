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

namespace palm_tracer {
/// <summary>Classe de calibration de la spline pour l'ajustement.</summary>
class Spline {
public:
	size_t roi_{},
	       sizeX_{},
	       sizeY_{},
	       sizeZ_{},
	       sizeXY_{},
	       sizeXYZ_{},
	       z0_{},
	       dz_{};

	Spline()  = default;
	~Spline() = default;

	void Calibrate(const uint16_t* image, const double* beads, size_t nBeads, size_t height, size_t width, size_t planes, size_t roi, size_t distance, size_t dz, size_t smoothz);

	void Load(const std::string& filename);

	void Save(const std::string& filename);


protected:
};
} // namespace palm_tracer
