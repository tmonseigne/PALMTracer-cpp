#include "TestUtils.hpp"
#include "signal_processing/PSF.hpp"

#include <array>
#include <numbers>

namespace palm_tracer::signal_processing {
namespace {
/// <summary>Position relative permettant de vérifier les valeurs et gradients de la PSF.</summary>
struct PSFCase {
	const char* name; ///< Identifiant du scénario.
	double x;         ///< Coordonnee relative horizontale.
	double y;         ///< Coordonnee relative verticale.
	double theta;     ///< Angle en radians.
};

/// <summary>Affiche le nom du scénario plutot que les octets de sa structure dans les rapports GoogleTest.</summary>
void PrintTo(const PSFCase& value, std::ostream* stream) { *stream << value.name; }

/// <summary>Scénarios analytiques de PSF, indépendants des images de microscopie.</summary>
class PSFTests : public testing::TestWithParam<PSFCase> { };
} // namespace

#pragma region PSF
/// <summary>La gaussienne isotrope suit sa valeur analytique et est paire autour de son centre.</summary>
TEST_P(PSFTests, IsotropicValue)
{
	const auto& [name, x, y, theta] = GetParam();
	constexpr double Offset         = 7, Amplitude = 120, Sigma = 1.3;
	const double Expected           = Offset + (Amplitude * std::exp(-((x * x) + (y * y)) / (2 * Sigma * Sigma)));
	EXPECT_NEAR(PSF(Offset, Amplitude, x, y, Sigma), Expected, 1e-12);
	EXPECT_DOUBLE_EQ(PSF(Offset, Amplitude, x, y, Sigma), PSF(Offset, Amplitude, -x, -y, Sigma));
}

/// <summary>Des largeurs égales à sqrt(2) sigma retrouvent le modèle isotrope pour tout angle.</summary>
TEST_P(PSFTests, CircularAnisotropicMatchesIsotropic)
{
	const auto& [name, x, y, theta] = GetParam();
	constexpr double Sigma          = 1.3, Width = Sigma * std::numbers::sqrt2;
	EXPECT_NEAR(PSF(7, 120, x, y, Width, Width, theta), PSF(7, 120, x, y, Sigma), 1e-12);
}
#pragma endregion PSF

#pragma region Deriv
/// <summary>Compare les cinq dérivées isotropes à des différences finies centrées, y compris les signes de x0 et y0.</summary>
TEST_P(PSFTests, IsotropicGradientMatchesFiniteDifferences)
{
	const auto& [name, x, y, theta] = GetParam();
	constexpr std::array<double, 5> Params{7, 120, 0, 0, 1.3};
	std::vector<double> gradient(5, -12345);
	Deriv(gradient, Params[1], x, y, Params[4], true);
	const auto Evaluate   = [x, y](const auto& p) { return PSF(p[0], p[1], x - p[2], y - p[3], p[4]); };
	constexpr double Step = 1e-5;
	for (size_t i = 0; i < Params.size(); ++i) {
		auto plus = Params, minus = Params;
		plus[i] += Step;
		minus[i] -= Step;
		EXPECT_NEAR(gradient[i], (Evaluate(plus) - Evaluate(minus)) / (2 * Step), 1e-7) << "Parametre " << i;
	}
	std::vector<double> fixed(5, -12345);
	Deriv(fixed, Params[1], x, y, Params[4], false);
	EXPECT_DOUBLE_EQ(fixed.back(), -12345);
	for (size_t i = 0; i < 4; ++i) { EXPECT_DOUBLE_EQ(fixed[i], gradient[i]); }
}

/// <summary>Compare les sept dérivées du modèle orienté à une évaluation numérique indépendante.</summary>
TEST_P(PSFTests, AnisotropicGradientMatchesFiniteDifferences)
{
	const auto& [name, x, y, theta] = GetParam();
	const std::array<double, 7> Params{7, 120, 0, 0, 1.3, 2.1, theta};
	std::vector<double> gradient(7, -12345);
	Deriv(gradient, Params[1], x, y, Params[4], Params[5], Params[6], true);
	const auto Evaluate   = [x, y](const auto& p) { return PSF(p[0], p[1], x - p[2], y - p[3], p[4], p[5], p[6]); };
	constexpr double Step = 1e-5;
	for (size_t i = 0; i < Params.size(); ++i) {
		auto plus = Params, minus = Params;
		plus[i] += Step;
		minus[i] -= Step;
		EXPECT_NEAR(gradient[i], (Evaluate(plus) - Evaluate(minus)) / (2 * Step), 1e-7) << "Parametre " << i;
	}
	std::vector<double> fixed(7, -12345);
	Deriv(fixed, Params[1], x, y, Params[4], Params[5], Params[6], false);
	EXPECT_DOUBLE_EQ(fixed.back(), -12345);
	for (size_t i = 0; i < 6; ++i) { EXPECT_DOUBLE_EQ(fixed[i], gradient[i]); }
}

INSTANTIATE_TEST_SUITE_P(Positions, PSFTests,
                         testing::Values(PSFCase{"Center", 0, 0, 0}, PSFCase{"XAxis", 0.7, 0, 0}, PSFCase{"YAxis", 0, -0.8, 0},
	                         PSFCase{"Diagonal", 0.7, -0.8, 0}, PSFCase{"OneDegree", 0.7, -0.8, std::numbers::pi / 180},
	                         PSFCase{"Rotated", 0.7, -0.8, 0.4},
	                         PSFCase{"QuarterTurn", -0.7, 0.8, std::numbers::pi / 2}),
                         [](const testing::TestParamInfo<PSFCase>& info) { return info.param.name; });
#pragma endregion Deriv
} // namespace palm_tracer::signal_processing
