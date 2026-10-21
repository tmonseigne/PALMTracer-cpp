#include "Palm.hpp"
#include "TestUtils.hpp"
#include "TiffManagement.hpp"

#include <array>
#include <limits>
#include <numbers>
#include <omp.h>
#include <tuple>

namespace palm_tracer {
namespace {
/// <summary>Image et dimensions conservées ensemble pour les scénarios de localisation.</summary>
struct Stack {
	std::vector<uint16_t> pixels; ///< Pixels contigus, plan par plan.
	Vec3I dim{};                  ///< Dimensions de la pile.

	/// <summary>Charge une pile de référence sans modifier le fichier source.</summary>
	explicit Stack(const std::string& name) { LoadTiff3D((InputPath / name).string().c_str(), pixels, dim); }
};

/// <summary>Exécute la localisation et contrôle le contrat du tampon de sortie.</summary>
Table Localize(const uint16_t* pixels, const Vec3I& dim, const size_t fit, const bool watershed)
{
	constexpr std::array Params{7.0, 1.0, 2.0, 0.0};
	// Un emplacement par pixel couvre largement la densité des piles de référence.
	std::vector<double> output(((dim.GetPlaneSize() / 5) + 1) * dim.z * localization_columns::Count, -12345.0);
	const Palm Processor(pixels, output.data(), output.size(), dim, 103.6, watershed ? 0.0 : 10.0, fit, Params.data());
	const int Count = Processor.Processing();
	EXPECT_GE(Count, 0);
	if (Count < 0 || static_cast<size_t>(Count) > output.size()) {
		ADD_FAILURE() << "Taille de sortie invalide";
		return {};
	}
	const auto Size = static_cast<size_t>(Count);
	EXPECT_EQ(Size % localization_columns::Count, 0);
	if (Size < output.size()) { EXPECT_DOUBLE_EQ(output[Size], -12345.0); }
	Table result;
	for (size_t i = 0; i + localization_columns::Count <= Size; i += localization_columns::Count) {
		result.emplace_back(output.begin() + static_cast<ptrdiff_t>(i), output.begin() + static_cast<ptrdiff_t>(i + localization_columns::Count));
		// Les CSV historiques stockent l'angle en degrés, ramene a [-90, 90[.
		auto& row = result.back();
		if (fit == 4 && row[8] > 0 && row[9] > 0) {
			constexpr double Pi = std::numbers::pi;
			double angle        = std::fmod(row[10] + (Pi / 2), Pi);
			if (angle < 0) { angle += Pi; }
			row[10] = (angle - (Pi / 2)) * 180 / Pi;
		}
	}
	return result;
}

/// <summary>Apparie chaque détection à une seule référence du même plan et canal, puis compare toutes les valeurs.</summary>
/// <remarks>Reprend les tolérances de compare_points : atol configurable et rtol de numpy.isclose à 1e-5.</remarks>
void ExpectLocalizationsNear(const Table& actual, const Table& expected, const double tolerance = 1e-5)
{
	ASSERT_EQ(actual.size(), expected.size());
	std::vector used(expected.size(), false);
	for (const auto& Row : actual) {
		SCOPED_TRACE("Plan " + std::to_string(Row[1]) + ", detection " + std::to_string(Row[0]));
		size_t nearest  = expected.size();
		double distance = std::numeric_limits<double>::max();
		for (size_t j = 0; j < expected.size(); ++j) {
			if (used[j] || Row[1] != expected[j][1] || Row[3] != expected[j][3]) { continue; }
			const double Dx = Row[4] - expected[j][4], Dy = Row[5] - expected[j][5], Dz = Row[6] - expected[j][6];
			const double D  = (Dx * Dx) + (Dy * Dy) + (Dz * Dz);
			if (D < distance) {
				nearest  = j;
				distance = D;
			}
		}
		ASSERT_LT(nearest, expected.size());
		used[nearest] = true;
		ExpectTableNear({Row}, {expected[nearest]}, tolerance, 1e-5);
	}
}

using LocalizationParams = std::tuple<size_t, bool>;

/// <summary>Variantes indépendantes des cinq modes disposant de références CSV et du watershed.</summary>
class PalmLocalizationTests : public testing::TestWithParam<LocalizationParams> { };

/// <summary>Restaure le budget OpenMP après chaque scénario, même en cas d'échec d'assertion.</summary>
struct ThreadBudget {
	int previous = omp_get_max_threads(); ///< Budget de l'appelant à restaurer.
	explicit ThreadBudget(const int threads) { omp_set_num_threads(threads); }
	~ThreadBudget() { omp_set_num_threads(previous); }
};

/// <summary>Comparaison du pipeline parallèle avec l'analyse indépendante des frames.</summary>
class PalmParallelTests : public testing::TestWithParam<std::tuple<size_t, int>> { };

/// <summary>Capacité réservée par plan, y compris avec une frame initiale sans détection.</summary>
class PalmCapacityTests : public testing::TestWithParam<bool> { };

/// <summary>Construit le suffixe historique commun aux fichiers Python et C++.</summary>
std::string LocalizationSuffix(const size_t fit, const bool watershed)
{
	return "103.6_" + std::string(watershed ? "True_" : "False_") + std::to_string(fit) + "_1.0_0.0_7.csv";
}
} // namespace

#pragma region Localisation
/// <summary>Une petite sortie réserve un emplacement à chaque plan et ne déborde pas.</summary>
TEST(PalmTests, BoundedOutput)
{
	const Stack Input("stack.tif");
	constexpr std::array Params{7.0, 1.0, 2.0, 0.0};
	const size_t Capacity = Input.dim.z * localization_columns::Count;
	std::vector<double> output(Capacity + 2, -12345);
	const Palm Processor(Input.pixels.data(), output.data() + 1, Capacity, Input.dim, 103.6, 0, 0, Params.data());
	EXPECT_EQ(Processor.Processing(), static_cast<int>(Capacity));
	EXPECT_DOUBLE_EQ(output.front(), -12345);
	EXPECT_DOUBLE_EQ(output.back(), -12345);
	for (size_t i = 0; i < Input.dim.z; ++i) {
		EXPECT_DOUBLE_EQ(output[1 + (i * localization_columns::Count)], static_cast<double>(i + 1));
		EXPECT_DOUBLE_EQ(output[1 + (i * localization_columns::Count) + 1], static_cast<double>(i + 1));
		EXPECT_DOUBLE_EQ(output[1 + (i * localization_columns::Count) + 2], 1);
		EXPECT_GT(output[1 + (i * localization_columns::Count) + 16], 0);
	}
}

/// <summary>Une frame vide ne transfère pas son quota aux autres ; une fraction de ligne reste inutilisée.</summary>
TEST_P(PalmCapacityTests, RespectsPerFrameCapacity)
{
	const ThreadBudget Budget(2);
	const Stack Input("stack.tif");
	const size_t PlaneSize = Input.dim.GetPlaneSize();
	const Vec3I SingleDim{.x = Input.dim.x, .y = Input.dim.y, .z = 1};
	const auto Expected = Localize(Input.pixels.data(), SingleDim, 0, true);
	ASSERT_GT(Expected.size(), 1);
	const Vec3I Dim{.x = Input.dim.x, .y = Input.dim.y, .z = GetParam() ? size_t{3} : size_t{2}};
	std::vector<uint16_t> pixels(Dim.GetSize(), 0);
	const size_t First = GetParam() ? 1 : 0;
	for (size_t i = First; i < Dim.z; ++i) { std::copy_n(Input.pixels.begin(), PlaneSize, pixels.begin() + static_cast<ptrdiff_t>(i * PlaneSize)); }
	constexpr std::array Params{7.0, 1.0, 2.0, 0.0};
	constexpr size_t PerPlane = 2;
	const size_t Rows         = (Dim.z - First) * PerPlane;
	const size_t Capacity     = (Dim.z * PerPlane * localization_columns::Count) + 1;
	std::vector<double> output(Capacity + 2, -12345);
	const Palm Processor(pixels.data(), output.data() + 1, Capacity, Dim, 103.6, 0, 0, Params.data());
	EXPECT_EQ(Processor.Processing(), static_cast<int>(Rows * localization_columns::Count));
	EXPECT_DOUBLE_EQ(output.front(), -12345);
	EXPECT_DOUBLE_EQ(output[1 + (Rows * localization_columns::Count)], -12345);
	EXPECT_DOUBLE_EQ(output[Capacity], -12345);
	EXPECT_DOUBLE_EQ(output.back(), -12345);
	for (size_t i = 0; i < Rows; ++i) {
		const size_t Offset = 1 + (i * localization_columns::Count);
		EXPECT_DOUBLE_EQ(output[Offset], static_cast<double>(i + 1));
		EXPECT_DOUBLE_EQ(output[Offset + 1], static_cast<double>(First + (i / PerPlane) + 1));
		EXPECT_DOUBLE_EQ(output[Offset + 2], static_cast<double>((i % PerPlane) + 1));
		EXPECT_GT(output[Offset + 16], 0);
	}
}

INSTANTIATE_TEST_SUITE_P(, PalmCapacityTests, testing::Bool(),
                           [](const testing::TestParamInfo<bool>& info) { return info.param ? "LeadingEmptyFrame" : "DenseFirstFrame"; });

/// <summary>Les frames vides et le budget de threads ne changent ni les fits ni l'ordre des sorties.</summary>
TEST_P(PalmParallelTests, MatchesIndividualFrames)
{
	const auto [fit, threads] = GetParam();
	const ThreadBudget Budget(threads);
	const Stack Input("stack.tif");
	const Vec3I Dim{.x = Input.dim.x, .y = Input.dim.y, .z = 5};
	const size_t PlaneSize = Dim.GetPlaneSize();
	std::vector<uint16_t> pixels(Dim.GetSize(), 0);
	// Des frames vides avant, entre et après les deux frames utiles vérifient l'assemblage des sorties.
	for (const size_t Plane : {size_t{1}, size_t{3}}) {
		std::copy_n(Input.pixels.begin(), PlaneSize, pixels.begin() + static_cast<ptrdiff_t>(Plane * PlaneSize));
	}
	const Vec3I SingleDim{.x = Dim.x, .y = Dim.y, .z = 1};
	const auto Single = Localize(Input.pixels.data(), SingleDim, fit, true);
	ASSERT_FALSE(Single.empty());
	const auto Actual = Localize(pixels.data(), Dim, fit, true);
	ASSERT_EQ(Actual.size(), 2 * Single.size());
	for (size_t i = 0; i < Actual.size(); ++i) {
		EXPECT_DOUBLE_EQ(Actual[i][0], static_cast<double>(i + 1));
		EXPECT_DOUBLE_EQ(Actual[i][1], i < Single.size() ? 2 : 4);
		EXPECT_DOUBLE_EQ(Actual[i][2], static_cast<double>((i % Single.size()) + 1));
		ExpectTableNear({Table::value_type(Actual[i].begin() + 3, Actual[i].end())},
		                {Table::value_type(Single[i % Single.size()].begin() + 3, Single[i % Single.size()].end())}, 1e-12);
	}
}

INSTANTIATE_TEST_SUITE_P(, PalmParallelTests, testing::Combine(testing::Values(size_t{0}, size_t{4}), testing::Values(1, 2)),
                           [](const testing::TestParamInfo<std::tuple<size_t, int>>& info) {
                           return std::string(std::get<0>(info.param) == 0 ? "NoFit" : "GaussianTheta") + "Threads" + std::to_string(std::get<1>(info.param));
                           });

/// <summary>Plusieurs frames différentes et vides réutilisent chaque espace sans transférer de détection ni de fit.</summary>
TEST_P(PalmParallelTests, ReusesWorkspacesAcrossDifferentFrames)
{
	const auto [fit, threads] = GetParam();
	const ThreadBudget Budget(threads);
	const Stack Input("stack.tif");
	const Vec3I Dim{.x = Input.dim.x, .y = Input.dim.y, .z = 9};
	const Vec3I SingleDim{.x = Dim.x, .y = Dim.y, .z = 1};
	const size_t PlaneSize = Dim.GetPlaneSize();
	std::vector<uint16_t> pixels(Dim.GetSize(), 0);
	Table expected;
	for (size_t plane = 0; plane < Dim.z; plane += 2) {
		const auto* image = Input.pixels.data() + ((plane % Input.dim.z) * PlaneSize);
		std::copy_n(image, PlaneSize, pixels.data() + (plane * PlaneSize));
		auto single = Localize(image, SingleDim, fit, true);
		for (auto& row : single) {
			row[0] = static_cast<double>(expected.size() + 1);
			row[1] = static_cast<double>(plane + 1);
			expected.push_back(row);
		}
	}
	const auto Actual = Localize(pixels.data(), Dim, fit, true);
	ASSERT_EQ(Actual.size(), expected.size());
	for (size_t row = 0; row < Actual.size(); ++row) { ExpectTableNear({Actual[row]}, {expected[row]}, 1e-12); }
}

/// <summary>Une capacité nulle ou une pile sans plan ne lit ni n'ecrit de tampon.</summary>
TEST(PalmTests, EmptyOutput)
{
	constexpr std::array Params{7.0, 1.0, 2.0, 0.0};
	const Palm NoCapacity(nullptr, nullptr, 0, {.x = 16, .y = 16, .z = 2}, 0, 0, 0, Params.data());
	EXPECT_EQ(NoCapacity.Processing(), 0);
	const Palm NoPlanes(nullptr, nullptr, 0, {.x = 16, .y = 16, .z = 0}, 0, 0, 0, Params.data());
	EXPECT_EQ(NoPlanes.Processing(), 0);
}

/// <summary>Une pile noire ne produit aucune détection et laisse le tampon intact.</summary>
TEST(PalmTests, BlackStack)
{
	constexpr Vec3I Dim{.x = 16, .y = 12, .z = 3};
	const std::vector<uint16_t> Pixels(Dim.GetSize(), 0);
	EXPECT_TRUE(Localize(Pixels.data(), Dim, 0, true).empty());
}

/// <summary>Compare tous les champs de localisation aux valeurs attendues sur la pile.</summary>
TEST_P(PalmLocalizationTests, Stack)
{
	const auto [fit, watershed] = GetParam();
	const Stack Input("stack.tif");
	const auto Expected = LoadCsv(InputPath / "ref" / ("stack-localizations-" + LocalizationSuffix(fit, watershed)));
	ExpectLocalizationsNear(Localize(Input.pixels.data(), Input.dim, fit, watershed), Expected, 1e-3);
}

INSTANTIATE_TEST_SUITE_P(, PalmLocalizationTests,
                           testing::Combine(testing::Values(size_t{0}, size_t{1}, size_t{2}, size_t{3}, size_t{4}), testing::Bool()),
                           ([](const testing::TestParamInfo<LocalizationParams>& info) {
	                           constexpr std::array Names{"NoFit", "GaussianXY", "GaussianSigma", "GaussianSigmaXY", "GaussianTheta"};
	                           return std::string(Names[std::get<0>(info.param)]) + (std::get<1>(info.param) ? "WithWatershed" : "WithoutWatershed");
	                           }));

/// <summary>Vérifie les quadrants éclairés et toutes les détections de la pile masquee.</summary>
TEST(PalmTests, Quadrants)
{
	const Stack Input("stack_quadrant.tif");
	const auto Actual = Localize(Input.pixels.data(), Input.dim, 4, true);
	ExpectLocalizationsNear(Actual, LoadCsv(InputPath / "ref/stack_quadrant-localizations-103.6_True_4_1.0_0.0_7.csv"));
	for (const auto& Row : Actual) {
		const auto Plane = static_cast<size_t>(Row[1]);
		const bool Top   = Plane == 3 || Plane == 4 || Plane == 7 || Plane == 8;
		const bool Left  = Plane == 1 || Plane == 4 || Plane == 5 || Plane == 8 || Plane == 9;
		EXPECT_TRUE(Top ? Row[5] <= 128 : Row[5] >= 128);
		EXPECT_TRUE(Left ? Row[4] <= 128 : Row[4] >= 128);
	}
}
#pragma endregion Localisation

#pragma region Seuil automatique
namespace {
/// <summary>Un cas indépendant par plan du calcul de seuil automatique.</summary>
class AutoThresholdTests : public testing::TestWithParam<size_t> { };
} // namespace

/// <summary>Compare le seuil de chaque plan à sa valeur attendue.</summary>
TEST_P(AutoThresholdTests, Value)
{
	constexpr std::array Expected{
		63.651749, 66.142910, 63.192296, 64.368855, 63.954150,
		63.400043, 66.521994, 63.373237, 62.613244, 63.866017
	};
	constexpr std::array Params{7.0};
	const Stack Input("stack.tif");
	ASSERT_LT(GetParam(), Input.dim.z);
	const Vec3I Dim{.x = Input.dim.x, .y = Input.dim.y, .z = 1};
	std::vector<double> output(((Dim.GetPlaneSize() / 5) + 1) * localization_columns::Count);
	Palm processor(Input.pixels.data() + (GetParam() * Dim.GetPlaneSize()), output.data(), output.size(), Dim, 0, 10, 0, Params.data());
	EXPECT_NEAR(processor.AutoThreshold(), Expected[GetParam()], 1e-5);
}

INSTANTIATE_TEST_SUITE_P(, AutoThresholdTests, testing::Range(size_t{0}, size_t{10}),
                           [](const testing::TestParamInfo<size_t>& info) { return "Plane" + std::to_string(info.param + 1); });
#pragma endregion Seuil automatique
} // namespace palm_tracer
