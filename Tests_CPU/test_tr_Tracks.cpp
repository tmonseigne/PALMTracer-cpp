#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <numeric>

#include "Structs.hpp"
#include "tracking/Tracks.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TEST(TracksTests, BlinkingReconnection)
{
	const std::vector<double> InitialTracks = {
		1, 1, 1, 1, 1, 1, 1, 1, 1, 99, 2, 1, 1, 1, 2, 2,
		2, 1, 3, 6, 1, 1, 3, 3, 2, 2, 4, 6, 1, 1, 4, 4,
		3, 3, 5, 4, 1, 1, 5, 5, 3, 4, 6, 4, 1, 1, 6, 6,
		4, 3, 7, 5, 1, 1, 7, 7, 4, 4, 8, 5, 1, 1, 8, 8,
		5, 5, 9, 4, 1, 1, 9, 9, 5, 6, 10, 4, 1, 1, 10, 10,
		6, 11, 11, 1, 1, 1, 11, 11, 6, 12, 12, 1, 1, 1, 12, 12,
		7, 16, 13, 5, 1, 1, 13, 13, 7, 17, 14, 5, 1, 1, 14, 14,
		8, 31, 15, 1, 1, 1, 15, 15, 8, 32, 16, 1, 1, 1, 16, 16,
		9, 36, 17, 9, 1, 1, 17, 17, 9, 37, 18, 9, 1, 1, 18, 18
	};

	const std::array Ref = {
		std::vector<double>{
			1, 1, 1, 1, 1, 1, 1, 1, 1, 99, 2, 1, 1, 1, 2, 2,
			2, 1, 3, 6, 1, 1, 3, 3, 2, 2, 4, 6, 1, 1, 4, 4,
			2, 3, 7, 5, 1, 1, 7, 7, 2, 4, 8, 5, 1, 1, 8, 8,
			2, 5, 9, 4, 1, 1, 9, 9, 2, 6, 10, 4, 1, 1, 10, 10,
			3, 3, 5, 4, 1, 1, 5, 5, 3, 4, 6, 4, 1, 1, 6, 6,
			4, 11, 11, 1, 1, 1, 11, 11, 4, 12, 12, 1, 1, 1, 12, 12,
			5, 16, 13, 5, 1, 1, 13, 13, 5, 17, 14, 5, 1, 1, 14, 14,
			6, 31, 15, 1, 1, 1, 15, 15, 6, 32, 16, 1, 1, 1, 16, 16,
			7, 36, 17, 9, 1, 1, 17, 17, 7, 37, 18, 9, 1, 1, 18, 18
		},
		std::vector<double>{
			1, 1, 1, 1, 1, 1, 1, 1, 1, 99, 2, 1, 1, 1, 2, 2,
			2, 1, 3, 6, 1, 1, 3, 3, 2, 2, 4, 6, 1, 1, 4, 4,
			2, 3, 7, 5, 1, 1, 7, 7, 2, 4, 8, 5, 1, 1, 8, 8,
			2, 5, 9, 4, 1, 1, 9, 9, 2, 6, 10, 4, 1, 1, 10, 10,
			3, 3, 5, 4, 1, 1, 5, 5, 3, 4, 6, 4, 1, 1, 6, 6,
			4, 11, 11, 1, 1, 1, 11, 11, 4, 12, 12, 1, 1, 1, 12, 12,
			4, 16, 13, 5, 1, 1, 13, 13, 4, 17, 14, 5, 1, 1, 14, 14,
			5, 31, 15, 1, 1, 1, 15, 15, 5, 32, 16, 1, 1, 1, 16, 16,
			6, 36, 17, 9, 1, 1, 17, 17, 6, 37, 18, 9, 1, 1, 18, 18
		},
		std::vector<double>{
			1, 1, 1, 1, 1, 1, 1, 1, 1, 99, 2, 1, 1, 1, 2, 2,
			2, 1, 3, 6, 1, 1, 3, 3, 2, 2, 4, 6, 1, 1, 4, 4,
			2, 3, 7, 5, 1, 1, 7, 7, 2, 4, 8, 5, 1, 1, 8, 8,
			2, 5, 9, 4, 1, 1, 9, 9, 2, 6, 10, 4, 1, 1, 10, 10,
			3, 3, 5, 4, 1, 1, 5, 5, 3, 4, 6, 4, 1, 1, 6, 6,
			4, 11, 11, 1, 1, 1, 11, 11, 4, 12, 12, 1, 1, 1, 12, 12,
			4, 16, 13, 5, 1, 1, 13, 13, 4, 17, 14, 5, 1, 1, 14, 14,
			5, 31, 15, 1, 1, 1, 15, 15, 5, 32, 16, 1, 1, 1, 16, 16,
			5, 36, 17, 9, 1, 1, 17, 17, 5, 37, 18, 9, 1, 1, 18, 18
		}
	};

	const size_t Size            = InitialTracks.size();
	constexpr size_t MaxDuration = 4;
	constexpr double MaxSpeed    = 2;
	Tracks tracks;
	tracks.LoadTracking(InitialTracks.data(), 18, 1);
	std::cout << " Initial Tracks : \n" << tracks << "\n";

	std::vector<double> output(Size, 0);
	for (size_t mode = 0; mode < 3; ++mode) {
		tracks.BlinkingReconnection(mode, MaxDuration, MaxSpeed);
		std::cout << " Reconnected Tracks (" << mode << ", " << MaxDuration << ", " << MaxSpeed << ") : \n" << tracks << "\n";
		tracks.FillOutput(output.data());
		for (size_t i = 0; i < Size; ++i) { EXPECT_DOUBLE_EQ(output[i], Ref[mode][i]); }
	}
}

//-------------------------------------------------------------------------------------------------
TEST(TracksTests, MSD)
{
	const std::vector<double> InitialTracks = {
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 5, 2, 2, 1, 1, 2, 2,

		2, 1, 3, 1, 1, 1, 3, 3,
		2, 2, 4, 2, 2, 1, 4, 4,
		2, 3, 5, 3, 0, 3, 7, 7,
		2, 4, 6, 7, 4, 5, 8, 8,

		3, 1, 7, 1, 1, 1, 11, 11,
		3, 2, 8, 3, 1, 1, 12, 12,
		3, 6, 9, 8, 1, 1, 13, 13,
		3, 7, 10, 2, 1, 1, 14, 14,
	};
	const std::vector<std::vector<double>> Ref = {
		{-1, -1, -1, 1},
		{15.6666667, 27, 61},
		{20, -1, -1, 25, 25, 1},
	};
	Tracks tracks;
	tracks.LoadTracking(InitialTracks.data(), 10, 1);
	const auto Res = tracks.ComputeAllMSD();
	for (size_t i = 0; i < Ref.size(); ++i) { for (size_t j = 0; j < Ref.size(); ++j) { EXPECT_NEAR(Res[i][j], Ref[i][j], 1e-5); } }

	// Print
	for (const auto& R : Res) {
		std::cout << " | ";
		for (const auto& C : R) { std::cout << C << " | "; }
		std::cout << "\n";
	}
}

//-------------------------------------------------------------------------------------------------
TEST(TracksTests, InstantDiffusion)
{
	const std::vector<double> InitialTracks = {
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 2, 2, 2, 2, 2, 2, 2,
		1, 3, 3, 3, 3, 3, 3, 3,
		1, 4, 4, 4, 4, 4, 4, 4,
	};
	const std::vector<std::vector<double>> Ref = {{1.5, 1.5}};
	Tracks tracks;
	tracks.LoadTracking(InitialTracks.data(), 4, 1);
	const auto Res = tracks.ComputeAllInstantDiffusion(2, true);
	for (size_t i = 0; i < Ref.size(); ++i) { for (size_t j = 0; j < Ref[i].size(); ++j) { EXPECT_DOUBLE_EQ(Res[i][j], Ref[i][j]); } }

	// Print
	for (const auto& R : Res) {
		std::cout << " | ";
		for (const auto& C : R) { std::cout << C << " | "; }
		std::cout << "\n";
	}
}

TEST(TracksTests, ComputeAllFit)
{
	const std::vector<double> InitialTracks = {
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 2, 2, 2, 2, 2, 2, 2,
		1, 3, 3, 3, 3, 3, 3, 3,
		1, 4, 4, 4, 4, 4, 4, 4,
	};
	const std::vector<std::vector<std::vector<double>>> Ref = {
		{{1.5, -6, 1, 12, -10, 0.979592}},                                                 // linear
		{{1.5, -6, 1, 2, 3, 1, 1.299038}},                                                 // power
		{{1.5, -6, 1, 29012085.510261, 2417672.688433, -4.217218, 0.979592, 4664.660863}}, // exponential
	};
	Tracks tracks;
	tracks.LoadTracking(InitialTracks.data(), 4, 1);
	const std::vector Modes  = {fit_mode::linear, fit_mode::power, fit_mode::exponential};
	const std::vector Params = {2.0};
	for (size_t i = 0; i < Modes.size(); ++i) {
		std::cout << "Mode " << static_cast<int>(Modes[i]) << " : \n";
		const auto Res = tracks.ComputeAllFitMetrics(true, Modes[i], Params);
		for (size_t j = 0; j < Ref[i].size(); ++j) { for (size_t k = 0; k < Ref[i][j].size(); ++k) { EXPECT_NEAR(Res[j][k], Ref[i][j][k], 1e-5); } }

		// Print
		for (const auto& R : Res) {
			std::cout << " | ";
			for (const auto& C : R) { std::cout << C << " | "; }
			std::cout << "\n";
		}
	}
}
} // namespace palm_tracer
