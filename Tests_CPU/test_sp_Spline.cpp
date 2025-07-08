#include "gtest/gtest.h"

#include <algorithm>

#include "Structs.hpp"
#include "TiffManagement.hpp"
#include "palm/Palm.hpp"
#include "signal_processing/Spline.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TEST(SplineTests, BaseCreation) { Spline spline; }

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, Calibrate)
{
	Spline spline;
	constexpr size_t ROI      = 29;
	constexpr size_t HALF_ROI = 14;

	std::vector<uint16_t> input;
	size_t height, width, planes;
	// 1. Charger l'image des billes
	LoadTiff3D("C:/Git/PALMTracer-cpp/Tests_CPU/Input/Beads.tiff", input, height, width, planes);
	const size_t PlaneSize = height * width;
	const size_t HalfPlane = planes / 2;
	std::cout << "Size(" << height << ", " << width << ", " << planes << "), Plane Size = " << PlaneSize << ", HalfPlane = " << HalfPlane;

	// Calcul du seuil sur l'image centrale (Z0) ROI de 25 pour prendre les billes complètes.
	auto PointsListSize = static_cast<size_t>(0.2 * static_cast<double>(height * width) * N_LOCALIZATIONS);
	std::vector<double> pointsList(PointsListSize);
	std::array<double, 4> FitParams = {ROI, 1.0, 2.0, 0.0};
	auto palm = Palm(input.data() + PlaneSize * HalfPlane, pointsList.data(), PointsListSize, height, width, 1, 0, 10, 0, FitParams.data());
	const double Threshold = palm.AutoThreshold();
	std::cout << ", Threshold = " << Threshold << "\n";

	// Calcul de la position des billes sur l'image centrale (Z0)
	size_t nBeads = palm.Processing();
	nBeads /= N_LOCALIZATIONS;
	std::cout << nBeads << " Beads found " << "\n";
	std::vector<Coord> beads;
	beads.reserve(nBeads);
	for (size_t i = 0; i < nBeads; ++i) {
		const size_t Idx = i * N_LOCALIZATIONS,
		             X   = static_cast<size_t>(std::round(pointsList[Idx + 4])),
		             Y   = static_cast<size_t>(std::round(pointsList[Idx + 5]));
		std::cout << "Bead " << i + 1 << " (" << X << ", " << Y << ")\n";
		if (HALF_ROI <= X && X < width - HALF_ROI && HALF_ROI <= Y && Y < height - HALF_ROI) { beads.emplace_back(Coord{.x = X, .y = Y}); }
	}
	std::cout << beads.size() << "\n";
	spline.Calibrate(input.data(), beads, height, width, planes, ROI, 50, 10, 1);
}
} // namespace palm_tracer
