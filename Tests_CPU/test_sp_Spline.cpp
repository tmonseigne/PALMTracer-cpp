#include "gtest/gtest.h"

#include <algorithm>

#include "Structs.hpp"
#include "TiffManagement.hpp"
#include "palm/Palm.hpp"
#include "signal_processing/Matrix.hpp"
#include "signal_processing/Spline.hpp"

namespace palm_tracer {
constexpr size_t ROI      = 29,
                 HALF_ROI = 14,
                 DISTANCE = 50,
                 DZ       = 10,
                 SMOOTH_Z = 1;

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, BaseCreation) { Spline spline; }


//-------------------------------------------------------------------------------------------------
TEST(SplineTests, setMembers)
{
	Spline spline;
	std::vector<uint16_t> input;
	Coord3D dim{};
	// Charger l'image des billes
	LoadTiff3D("C:/Git/PALMTracer-cpp/Tests_CPU/Input/Beads.tiff", input, dim);

	spline.setMembers(dim, ROI, DISTANCE, DZ, SMOOTH_Z);
	EXPECT_EQ(spline.smallXIndices_.size(), 13);           // Les sous-groupes seront inférieur à la ROI de base.
	EXPECT_EQ(spline.smallZIndices_.size(), DISTANCE + 1); // Prend un nombre impair de points

	spline.setMembers(dim, 8, 2 * dim.z, 5, 2);
	EXPECT_EQ(spline.smallXIndices_.size(), 7);     // Prend un nombre impair de points en roi et ne peut descnedre en dessous de la ROI (impair) de base
	EXPECT_EQ(spline.smallZIndices_.size(), dim.z); // On dépasse le max en Z donc on prend toutes les frames
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, extractBeads)
{
	Spline spline;
	std::vector<uint16_t> input;
	Coord3D dim{};
	// Charger l'image des billes
	LoadTiff3D("C:/Git/PALMTracer-cpp/Tests_CPU/Input/Beads.tiff", input, dim);

	const std::vector<Coord> Beads = {
		{.x = 25, .y = 19},
		{.x = 93, .y = 27},
		{.x = 122, .y = 61},
		{.x = 43, .y = 70},
		{.x = 44, .y = 82},
		{.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};

	std::vector<std::vector<double>> Ref = {
		{
			5649, 7067, 5629, 6618, 8244, 6662, 5981, 7379, 6055,
			5651, 6824, 5960, 7046, 8412, 6794, 6364, 8189, 6576,
			5433, 6535, 6006, 6416, 8187, 6565, 6056, 7679, 6704,
		},
		{
			5172, 6799, 6313, 5498, 7832, 6742, 4807, 6545, 6802,
			5443, 6705, 6205, 6229, 7863, 7687, 5633, 7464, 6935,
			5622, 6803, 6558, 6988, 7941, 8200, 5355, 7157, 7274,
		},
		{
			6267, 7407, 5861, 6604, 8559, 7510, 5780, 7121, 6292,
			5694, 6783, 6116, 6647, 8450, 6786, 5986, 6788, 6532,
			6383, 7393, 6380, 6965, 8189, 7150, 5536, 6763, 6658,
		},
		{
			5110, 6170, 5530, 5957, 7782, 7695, 5076, 7362, 7384,
			4826, 6413, 6168, 5378, 7752, 7723, 5052, 7791, 7494,
			4745, 6437, 6047, 5484, 7453, 7434, 5292, 7714, 7202,
		},
		{
			4474, 5975, 5819, 5118, 7618, 8008, 5646, 7652, 7189,
			4179, 6187, 6068, 6124, 7656, 7854, 5109, 7613, 7913,
			4349, 6056, 5916, 5278, 8270, 7618, 5212, 8199, 7725,
		},
		{
			6254, 6562, 4899, 7842, 7769, 5238, 7192, 7590, 5075,
			6173, 5778, 4300, 7036, 8591, 5715, 7265, 7623, 5472,
			6202, 6059, 4423, 7564, 7794, 5399, 6926, 7364, 5367,
		},
		{
			64, 194, 154, 83, 56, 58, 72, 58, 80,
			80, 210, 135, 113, 152, 96, 77, 133, 114,
			149, 104, 116, 152, 68, 123, 57, 117, 77,
		}
	};

	spline.setMembers(dim, 3, 3, DZ, SMOOTH_Z);
	const auto BeadsStack = spline.extractBeads(input.data(), Beads);

	// Sélection de 3 plans centraux uniquement pour comparaison
	for (size_t b = 0; b < BeadsStack.size(); ++b) {
		size_t idx = 0;
		for (size_t i = 0; i < spline.smallZIndices_.size(); ++i) {
			const size_t ZIdx = spline.smallZIndices_[i] * spline.beadImSize_;
			for (size_t j = 0; j < spline.beadImSize_; ++j, ++idx) {
				EXPECT_DOUBLE_EQ(Ref[b][idx], BeadsStack[b][ZIdx + j]);
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, Calibrate)
{
	Spline spline;

	std::vector<uint16_t> input;
	Coord3D dim{};
	// Charger l'image des billes
	LoadTiff3D("C:/Git/PALMTracer-cpp/Tests_CPU/Input/Beads.tiff", input, dim);
	const size_t PlaneSize = dim.GetPlaneSize();
	const size_t HalfPlane = dim.z / 2;
	std::cout << "Size " << dim << ", Plane Size = " << PlaneSize << ", HalfPlane = " << HalfPlane;

	// Calcul du seuil sur l'image centrale (Z0) ROI de 25 pour prendre les billes complètes.
	const auto PointsListSize = static_cast<size_t>(0.2 * static_cast<double>(PlaneSize) * N_LOCALIZATIONS);
	std::vector<double> pointsList(PointsListSize);
	constexpr std::array<double, 4> FitParams = {ROI, 1.0, 2.0, 0.0};
	auto palm = Palm(input.data() + PlaneSize * HalfPlane, pointsList.data(), PointsListSize, dim.y, dim.x, 1, 0, 10, 0, FitParams.data());
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
		if (HALF_ROI <= X && X < dim.x - HALF_ROI && HALF_ROI <= Y && Y < dim.y - HALF_ROI) { beads.emplace_back(Coord{.x = X, .y = Y}); }
	}
	std::cout << beads.size() << "\n";
	//spline.Calibrate(input.data(), beads, dim, ROI, 50, 10, 1);
}
} // namespace palm_tracer
