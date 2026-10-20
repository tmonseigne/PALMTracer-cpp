#include "TestUtils.hpp"
#include "fit/GaussianModel.hpp"
#include "signal_processing/PSF.hpp"

#include <array>
#include <stdexcept>
#include <tuple>

namespace palm_tracer::fit {
namespace {
/// <summary>Modèles isotropes/anisotropes avec leur dernier paramètre fixé/libre.</summary>
class GaussianModelTests : public testing::TestWithParam<std::tuple<bool, bool>> { };
} // namespace

/// <summary>Un état unique restitue valeurs et gradients sur plusieurs pixels sans état mutable résiduel.</summary>
TEST_P(GaussianModelTests, SharedValueAndGradientAcrossPixels)
{
	const auto [isotropic, variable] = GetParam();
	const GaussianModel Model        = isotropic
		                                   ? GaussianModel(7, 120, 1.3, variable)
		                                   : GaussianModel(7, 120, 1.3, 2.1, 0.4, variable);
	for (const auto& Position : std::array{std::array{0.0, 0.0}, std::array{0.7, -0.8}, std::array{-1.4, 0.2}, std::array{0.0, 0.0}}) {
		const double X = Position[0], Y = Position[1];
		std::vector<double> expected(8, -12345), actual(8, -12345);
		const double Value = isotropic ? PSF(7, 120, X, Y, 1.3) : PSF(7, 120, X, Y, 1.3, 2.1, 0.4);
		if (isotropic) { Deriv(expected, 120, X, Y, 1.3, variable); }
		else { Deriv(expected, 120, X, Y, 1.3, 2.1, 0.4, variable); }
		EXPECT_DOUBLE_EQ(Model.Evaluate(X, Y, actual), Value);
		EXPECT_DOUBLE_EQ(Model.Evaluate(X, Y), Value);
		for (size_t i = 0; i < expected.size(); ++i) { EXPECT_DOUBLE_EQ(actual[i], expected[i]) << "Composante " << i; }
	}
}

/// <summary>Une vue de gradient insuffisante est refusée avant toute écriture hors limites.</summary>
TEST_P(GaussianModelTests, RejectsShortGradient)
{
	const auto [isotropic, variable] = GetParam();
	const GaussianModel Model        = isotropic
		                                   ? GaussianModel(7, 120, 1.3, variable)
		                                   : GaussianModel(7, 120, 1.3, 2.1, 0.4, variable);
	const size_t Required = isotropic ? (variable ? 5 : 4) : (variable ? 7 : 6);
	std::vector<double> gradient(Required - 1, -12345);
	EXPECT_THROW(static_cast<void>(Model.Evaluate(0.7, -0.8, gradient)), std::invalid_argument);
	for (const double Value : gradient) { EXPECT_DOUBLE_EQ(Value, -12345); }
}

INSTANTIATE_TEST_SUITE_P(Models, GaussianModelTests, testing::Combine(testing::Bool(), testing::Bool()),
                         [](const testing::TestParamInfo<std::tuple<bool, bool>>& info) {
                         return std::string(std::get<0>(info.param) ? "Isotropic" : "Anisotropic") + (std::get<1>(info.param) ? "Free" : "Fixed");
                         });
} // namespace palm_tracer::fit
