#include "Palm.hpp"
#include "TestUtils.hpp"
#include "TiffManagement.hpp"

#include <array>
#include <limits>
#include <numbers>
#include <tuple>

namespace palm_tracer {
namespace {
/// <summary>Image et dimensions conservees ensemble pour les scenarios de localisation.</summary>
struct Stack {
	std::vector<uint16_t> pixels; ///< Pixels contigus, plan par plan.
	Vec3I dim{};                  ///< Dimensions de la pile.

	/// <summary>Charge une pile de reference sans modifier le fichier source.</summary>
	explicit Stack(const std::string& name) { LoadTiff3D((InputPath / name).string().c_str(), pixels, dim); }
};

/// <summary>Execute la localisation et controle le contrat du tampon de sortie.</summary>
Table Localize(const uint16_t* pixels, const Vec3I& dim, const size_t fit, const bool watershed)
{
	constexpr std::array Params{7.0, 1.0, 2.0, 0.0};
	// Un emplacement par pixel couvre largement la densite des piles de reference.
	std::vector<double> output((dim.GetPlaneSize() / 5 + 1) * dim.z * N_LOCALIZATIONS, -12345.0);
	const Palm Processor(pixels, output.data(), output.size(), dim, 103.6, watershed ? 0.0 : 10.0, fit, Params.data());
	const int Count = Processor.Processing();
	EXPECT_GE(Count, 0);
	if (Count < 0 || static_cast<size_t>(Count) > output.size()) {
		ADD_FAILURE() << "Taille de sortie invalide";
		return {};
	}
	const auto Size = static_cast<size_t>(Count);
	EXPECT_EQ(Size % N_LOCALIZATIONS, 0);
	if (Size < output.size()) { EXPECT_DOUBLE_EQ(output[Size], -12345.0); }
	Table result;
	for (size_t i = 0; i + N_LOCALIZATIONS <= Size; i += N_LOCALIZATIONS) {
		result.emplace_back(output.begin() + static_cast<ptrdiff_t>(i), output.begin() + static_cast<ptrdiff_t>(i + N_LOCALIZATIONS));
		// Les CSV historiques stockent l'angle en degres, ramene a [-90, 90[.
		auto& row = result.back();
		if (fit == 4 && row[8] > 0 && row[9] > 0) {
			constexpr double Pi = std::numbers::pi;
			double angle        = std::fmod(row[10] + Pi / 2, Pi);
			if (angle < 0) { angle += Pi; }
			row[10] = (angle - Pi / 2) * 180 / Pi;
		}
	}
	return result;
}

/// <summary>Apparie chaque detection a une seule reference du meme plan et canal, puis compare toutes les valeurs.</summary>
/// <remarks>Reprend les tolerances de compare_points : atol configurable et rtol de numpy.isclose a 1e-5.</remarks>
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
			const double D  = Dx * Dx + Dy * Dy + Dz * Dz;
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

/// <summary>Variantes independantes des cinq modes disposant de references CSV et du watershed.</summary>
class PalmLocalizationTests : public testing::TestWithParam<LocalizationParams> { };

/// <summary>Construit le suffixe historique commun aux fichiers Python et C++.</summary>
std::string LocalizationSuffix(const size_t fit, const bool watershed)
{
	return "103.6_" + std::string(watershed ? "True_" : "False_") + std::to_string(fit) + "_1.0_0.0_7.csv";
}
} // namespace

#pragma region Localisation
/// <summary>Une pile noire ne produit aucune detection et laisse le tampon intact.</summary>
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

/// <summary>Verifie les quadrants eclaires et toutes les detections de la pile masquee.</summary>
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
/// <summary>Un cas independant par plan du calcul de seuil automatique.</summary>
class AutoThresholdTests : public testing::TestWithParam<size_t> { };
}

/// <summary>Compare le seuil de chaque plan a sa valeur attendue.</summary>
TEST_P(AutoThresholdTests, Value)
{
	constexpr std::array Expected{
		63.639888, 65.789447, 63.192296, 64.375352, 63.954150,
		63.400043, 66.521994, 63.373237, 62.515444, 63.866017
	};
	constexpr std::array Params{7.0};
	const Stack Input("stack.tif");
	ASSERT_LT(GetParam(), Input.dim.z);
	const Vec3I Dim{.x = Input.dim.x, .y = Input.dim.y, .z = 1};
	std::vector<double> output((Dim.GetPlaneSize() / 5 + 1) * N_LOCALIZATIONS);
	Palm processor(Input.pixels.data() + GetParam() * Dim.GetPlaneSize(), output.data(), output.size(), Dim, 0, 10, 0, Params.data());
	EXPECT_NEAR(processor.AutoThreshold(), Expected[GetParam()], 1e-5);
}

INSTANTIATE_TEST_SUITE_P(, AutoThresholdTests, testing::Range(size_t{0}, size_t{10}),
                           [](const testing::TestParamInfo<size_t>& info) { return "Plane" + std::to_string(info.param + 1); });
#pragma endregion Seuil automatique
} // namespace palm_tracer
