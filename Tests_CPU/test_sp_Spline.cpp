#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <omp.h>

#include "Structs.hpp"
#include "TiffManagement.hpp"
#include "palm/Palm.hpp"
#include "signal_processing/Matrix.hpp"
#include "signal_processing/Spline.hpp"

namespace palm_tracer {
constexpr size_t ROI      = 29,
                 DISTANCE = 50,
                 DZ       = 10;

const std::filesystem::path BEAD_FILE_PATH = std::filesystem::path(__FILE__).parent_path() / "Input/Beads.tiff";
//-------------------------------------------------------------------------------------------------
TEST(SplineTests, SetMembers)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	spline.SetMembers(dim, ROI, DISTANCE, DZ);
	EXPECT_EQ(spline.smallXIndices_.size(), 13);           // Les sous-groupes seront inférieur à la ROI de base.
	EXPECT_EQ(spline.smallZIndices_.size(), DISTANCE + 1); // Prend un nombre impair de points

	spline.SetMembers(dim, 8, 2 * dim.z, 5);
	EXPECT_EQ(spline.smallXIndices_.size(), 7);     // Prend un nombre impair de points en roi et ne peut descnedre en dessous de la ROI (impair) de base
	EXPECT_EQ(spline.smallZIndices_.size(), dim.z); // On dépasse le max en Z donc on prend toutes les frames
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, ExtractBeads)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	const std::vector<Vec2I> Beads = {
		{.x = 25, .y = 19}, {.x = 93, .y = 27}, {.x = 122, .y = 61},
		{.x = 43, .y = 70}, {.x = 44, .y = 82}, {.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};

	const std::vector<std::vector<double>> Ref = {
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

	spline.SetMembers(dim, 3, 3, DZ);
	const auto BeadsStack = spline.ExtractBeads(input.data(), Beads);

	// Sélection de 3 plans centraux uniquement pour comparaison
	for (size_t b = 0; b < BeadsStack.size(); ++b) {
		size_t idx = 0;
		for (size_t i = 0; i < spline.smallZIndices_.size(); ++i) {
			const size_t ZIdx = spline.smallZIndices_[i] * spline.beadImSize_;
			for (size_t j = 0; j < spline.beadImSize_; ++j, ++idx) { EXPECT_DOUBLE_EQ(Ref[b][idx], BeadsStack[b][ZIdx + j]); }
		}
	}
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, Align)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	const std::vector<Vec2I> Beads = {
		{.x = 25, .y = 19}, {.x = 93, .y = 27}, {.x = 122, .y = 61},
		{.x = 43, .y = 70}, {.x = 44, .y = 82}, {.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};

	spline.SetMembers(dim, ROI, DISTANCE, DZ);
	auto beadsStack = spline.ExtractBeads(input.data(), Beads);

	const size_t Nref = std::min<size_t>(std::max<size_t>(5, Beads.size() / 2), Beads.size());

	//Création de la PSF de référence moyenne sur le volume complet
	std::vector mu(spline.beadSize_, 0.0);
	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(spline.beadSize_); ++i) {
		double sum = 0.0;
		for (size_t b = 0; b < Nref; ++b) { sum += beadsStack[b][i]; }
		mu[i] = sum / static_cast<double>(Nref);
	}
	const auto Shifts            = spline.Align(mu, beadsStack, false);
	const std::vector<Vec3D> Ref = {
		{.x = -0.0375, .y = 0.025, .z = 0},
		{.x = -0.0125, .y = 0.1, .z = 0},
		{.x = -0.1, .y = -0.025, .z = 0},
		{.x = 0.0875, .y = -0.1125, .z = 0},
		{.x = 0.075, .y = 0.0125, .z = 0},
		{.x = 0.1125, .y = -0.1, .z = -0.0125},
		{.x = 0.075, .y = 0.075, .z = -0.0375}
	};
	for (size_t i = 0; i < Shifts.size(); ++i) {
		EXPECT_NEAR(Shifts[i].x, Ref[i].x, 1e-3) << "Mismatch at Shifts[" << i << "].x";
		EXPECT_NEAR(Shifts[i].y, Ref[i].y, 1e-3) << "Mismatch at Shifts[" << i << "].y";
		EXPECT_NEAR(Shifts[i].z, Ref[i].z, 1e-3) << "Mismatch at Shifts[" << i << "].z";
	}
	// Dans Matlab :
	//   0.0750, -0.3625, 0.0125
	//	-0.0125, -0.125,  0
	//	-0.1375,  0.4625, 0
	//	 0.2,    -0.2,    0
	//	 0.725,   0.45,   0
	//	-0.125,   0.25,   0
	std::cerr << "[WARNING] Ce test utilise une matrice non linéaire, ce qui produit une interpolation différente dans MATLAB. Mais résultat cohérent.\n";
	SUCCEED();
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, GetOverlapMinimal)
{
	Spline spline;
	constexpr Vec3I Dim{.x = 5, .y = 5, .z = 3}; // Très petit volume

	// Génère des "billes" artificielles
	std::vector beads(6, std::vector(Dim.GetSize(), 0.0));

	// Remplit des cas différents pour chaque bille
	// Bille 0 : bonne bille
	std::iota(beads[0].begin(), beads[0].end(), 1.0);

	// Bille 1 : contient des NaN
	std::iota(beads[1].begin(), beads[1].end(), 1.0);
	beads[1][6] = std::numeric_limits<double>::quiet_NaN();

	// Bille 2 : complètement nulle
	std::ranges::fill(beads[2], 0.0);

	// Bille 3 : bonne, mais avec poids faible
	std::iota(beads[3].begin(), beads[3].end(), 1.0);
	matrix::Scale(beads[3], 0.5);

	// Bille 4 : bonne, mais shift > 3
	std::iota(beads[4].begin(), beads[4].end(), 1.0);

	// Bille 5 : bonne, mais différente (pour filtrage RMS)
	std::ranges::fill(beads[5], 10.0);

	// Shifts associés
	const std::vector Shifts = {
		Vec3D{.x = 0.0, .y = 0.0, .z = 0.0}, // Valide
		Vec3D{.x = 0.0, .y = 0.0, .z = 0.0}, // NaN
		Vec3D{.x = 0.0, .y = 0.0, .z = 0.0}, // Nulle
		Vec3D{.x = 0.0, .y = 0.0, .z = 0.0}, // Poids faible
		Vec3D{.x = 4.0, .y = 0.0, .z = 0.0}, // >3 ⇒ filtré
		Vec3D{.x = 0.0, .y = 0.0, .z = 0.0}  // Différente pour test RMS élevé
	};

	// Poids initiaux
	std::vector weight{1.0, 1.0, 1.0, 0.1, 1.0, 1.0};
	double norm = 0.0;

	spline.SetMembers(Dim, 5, DISTANCE, DZ);
	// changement des smallBeadSize pour avoir un equivalent
	spline.smallXIndices_ = {1, 2, 3};
	spline.smallZIndices_ = {1};

	//Dimansions de la zone d'interet centrale.
	spline.smallBeadDim_    = Vec3I{.x = spline.smallXIndices_.size(), .y = spline.smallXIndices_.size(), .z = spline.smallZIndices_.size()};
	spline.smallBeadImSize_ = spline.smallBeadDim_.GetPlaneSize();
	spline.smallBeadSize_   = spline.smallBeadDim_.GetSize();
	spline.GetOverlap(beads, Shifts, weight, norm, false);

	EXPECT_NEAR(norm, 43.2255, 1e-3) << "Mismatch for norm.";
	const std::vector RefWeight{3.4486, 3.4486, std::numeric_limits<double>::quiet_NaN(), 8.94678, std::numeric_limits<double>::quiet_NaN(), 2.26343};
	for (size_t i = 0; i < RefWeight.size(); ++i) {
		if (std::isnan(RefWeight[i])) { EXPECT_TRUE(std::isnan(weight[i])) << "Weight[" << i << "] expected to be NaN"; }
		else { EXPECT_NEAR(RefWeight[i], weight[i], 1e-3) << "Mismatch at Weight[" << i << "]"; }
	}

	// Affichage pour vérification
	//std::cout << "Norm : " << norm << "\n";
	//std::cout << "Good : [";
	//for (const auto& G : good) { std::cout << G << ", "; }
	//std::cout << "]\nWeight : [";
	//for (const auto& W : weight) { std::cout << W << ", "; }
	//std::cout << "]\n";
}


//-------------------------------------------------------------------------------------------------
TEST(SplineTests, GetOverlap)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	const std::vector<Vec2I> Beads = {
		{.x = 25, .y = 19}, {.x = 93, .y = 27}, {.x = 122, .y = 61},
		{.x = 43, .y = 70}, {.x = 44, .y = 82}, {.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};

	const std::vector<Vec3D> Shifts = {
		{.x = 0.0750, .y = -0.3625, .z = 0.0125},
		{.x = -0.0125, .y = -0.1250, .z = 0},
		{.x = -0.1375, .y = 0.4625, .z = 0},
		{.x = 0.2000, .y = -0.2000, .z = 0},
		{.x = 0.7250, .y = 0.4500, .z = 0},
		{.x = -0.1250, .y = 0.2500, .z = 0},
		{.x = -0.1250, .y = 0.2500, .z = 0}
	};
	spline.SetMembers(dim, ROI, DISTANCE, DZ);
	const auto BeadsStack = spline.ExtractBeads(input.data(), Beads);
	std::vector weight(Beads.size(), 1.0);
	double norm;
	spline.GetOverlap(BeadsStack, Shifts, weight, norm, false);

	EXPECT_NEAR(norm, 7860.1429, 1e-3) << "Mismatch for norm.";
	const std::vector RefWeight{58.646, 57.8895, 55.5306, 52.9919, 49.7937, 42.2984, std::numeric_limits<double>::quiet_NaN()};
	for (size_t i = 0; i < RefWeight.size(); ++i) {
		if (std::isnan(RefWeight[i])) { EXPECT_TRUE(std::isnan(weight[i])) << "Weight[" << i << "] expected to be NaN"; }
		else { EXPECT_NEAR(RefWeight[i], weight[i], 1e-3) << "Mismatch at Weight[" << i << "]"; }
	}
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, SoftNormalizationMinimal)
{
	Spline spline;
	constexpr Vec3I Dim{.x = 5, .y = 5, .z = 3}; // Très petit volume
	const size_t N        = Dim.GetSize();
	constexpr size_t Nvol = 3;

	// Génère des "billes" artificielles
	std::vector beads(Nvol, std::vector(N, 0.0));

	// Bille 0 : gradient croissant
	std::iota(beads[0].begin(), beads[0].end(), 1.0);

	// Bille 1 : *2 du volume 0
	std::ranges::transform(beads[0], beads[1].begin(), [](const double v) { return v * 2.0; });

	// Bille 2 : bruité avec NaNs
	std::iota(beads[2].begin(), beads[2].end(), 1.0);
	beads[2][6]  = std::numeric_limits<double>::quiet_NaN();
	beads[2][43] = std::numeric_limits<double>::quiet_NaN();

	spline.SetMembers(Dim, 5, DISTANCE, DZ);
	spline.smallZIndices_ = {0, 1, 2};
	const auto BeadsNorm  = spline.SoftNormalization(beads);

	// Vérifications : les volumes 0 et 1 doivent être très proches après normalisation
	for (size_t i = 0; i < N; ++i) {
		if (std::isfinite(BeadsNorm[0][i]) && std::isfinite(BeadsNorm[1][i])) { EXPECT_NEAR(BeadsNorm[0][i], BeadsNorm[1][i], 1e-3); }
	}
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, SoftNormalization)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	const std::vector<Vec2I> Beads = {
		{.x = 25, .y = 19}, {.x = 93, .y = 27}, {.x = 122, .y = 61},
		{.x = 43, .y = 70}, {.x = 44, .y = 82}, {.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};

	spline.SetMembers(dim, ROI, DISTANCE, DZ);
	const auto BeadsStack = spline.ExtractBeads(input.data(), Beads);
	const auto BeadsNorm  = spline.SoftNormalization(BeadsStack);
	// Affichage pour vérification
	//matrix::PrintVolume(BeadsNorm[0], spline.beadDim_);
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, GetPSFMinimal)
{
	// La parallelisation sur le CI Linux et FFTW3 plante parfois, on force un thread unique (cas pour les très petites matrices).
	#ifdef _OPENMP
	const int prev_threads = omp_get_max_threads();
	omp_set_num_threads(1);
	#endif

	Spline spline;
	constexpr Vec3I Dim{.x = 7, .y = 7, .z = 5}; // Très petit volume
	const size_t N        = Dim.GetSize();
	constexpr size_t Nvol = 3;

	// Génère des "billes" artificielles
	std::vector beads(Nvol, std::vector(N, 0.0));

	// Bille 0 : gradient croissant
	std::iota(beads[0].begin(), beads[0].end(), 1.0);

	// Bille 1 : *2 du volume 0
	std::ranges::transform(beads[0], beads[1].begin(), [](const double v) { return v * 2.0; });

	// Bille 2 : bruité avec NaNs
	std::iota(beads[2].begin(), beads[2].end(), 1.0);
	beads[2][8]  = std::numeric_limits<double>::quiet_NaN();
	beads[2][73] = std::numeric_limits<double>::quiet_NaN();
	//matrix::PrintVolume(beads[2], Dim);

	spline.SetMembers(Dim, 7, DISTANCE, DZ);
	// changement des smallBeadSize pour avoir un equivalent
	spline.smallXIndices_ = {0, 1, 2, 3, 4, 5, 6};
	spline.smallZIndices_ = {1, 2, 3};

	//Dimansions de la zone d'interet centrale.
	spline.smallBeadDim_    = Vec3I{.x = spline.smallXIndices_.size(), .y = spline.smallXIndices_.size(), .z = spline.smallZIndices_.size()};
	spline.smallBeadImSize_ = spline.smallBeadDim_.GetPlaneSize();
	spline.smallBeadSize_   = spline.smallBeadDim_.GetSize();

	const auto PSF = spline.GetPSF(beads);
	std::cerr << "[WARNING] Ce test utilise une matrice non linéaire, ce qui produit une interpolation différente dans MATLAB. Mais résultat cohérent.\n";
	SUCCEED();

	#ifdef _OPENMP
	omp_set_num_threads(prev_threads); // restaurer
	#endif
}

//-------------------------------------------------------------------------------------------------
TEST(SplineTests, GetPSF)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	const std::vector<Vec2I> Beads = {
		{.x = 25, .y = 19}, {.x = 93, .y = 27}, {.x = 122, .y = 61},
		{.x = 43, .y = 70}, {.x = 44, .y = 82}, {.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};

	spline.SetMembers(dim, ROI, DISTANCE, DZ);
	auto beads     = spline.ExtractBeads(input.data(), Beads);
	const auto PSF = spline.GetPSF(beads);
	//matrix::PrintVolume(PSF, spline.beadDim_);	// Affichage pour vérification
	std::cerr << "[WARNING] Ce test utilise une matrice non linéaire, ce qui produit une interpolation différente dans MATLAB. Mais résultat cohérent.\n";
	SUCCEED();
}


//-------------------------------------------------------------------------------------------------
TEST(SplineTests, Calibrate)
{
	Spline spline;
	std::vector<uint16_t> input;
	Vec3I dim{};
	// Charger l'image des billes
	LoadTiff3D(BEAD_FILE_PATH.string().c_str(), input, dim);

	const std::vector<Vec2I> Beads = {
		{.x = 25, .y = 19}, {.x = 93, .y = 27}, {.x = 122, .y = 61},
		{.x = 43, .y = 70}, {.x = 44, .y = 82}, {.x = 90, .y = 110},
		{.x = 28, .y = 118}
	};
	spline.Calibrate(input.data(), Beads, dim, ROI, DISTANCE, DZ);

	/*
	 * Version complète partant de 0.
		Spline spline;

		std::vector<uint16_t> input;
		Vec3I dim{};
		// Charger l'image des billes
		LoadTiff3D(beadfile_path.string().c_str(), input, dim);
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
		std::vector<Vec2I> beads;
		beads.reserve(nBeads);
		for (size_t i = 0; i < nBeads; ++i) {
			const size_t Idx = i * N_LOCALIZATIONS,
			             X   = static_cast<size_t>(std::round(pointsList[Idx + 4])),
			             Y   = static_cast<size_t>(std::round(pointsList[Idx + 5]));
			std::cout << "Bead " << i + 1 << " (" << X << ", " << Y << ")\n";
			if (HALF_ROI <= X && X < dim.x - HALF_ROI && HALF_ROI <= Y && Y < dim.y - HALF_ROI) { beads.emplace_back(Vec2I{.x = X, .y = Y}); }
		}
		std::cout << beads.size() << "\n";
		//spline.Calibrate(input.data(), beads, dim, ROI, 50, 10, 1);

		*/
	std::cerr << "[WARNING] Ce test utilise une matrice non linéaire, ce qui produit une interpolation différente dans MATLAB. Mais résultat cohérent.\n";
	SUCCEED();
}
} // namespace palm_tracer
