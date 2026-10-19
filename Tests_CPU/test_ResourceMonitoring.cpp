#include "ResourceMonitoring.hpp"
#include "TestUtils.hpp"

#include <tuple>

namespace palm_tracer::monitoring {
namespace {
using CpuCase = std::tuple<double, double, unsigned int, double>;

/// <summary>Scénarios déterministes de normalisation par la capacité CPU totale.</summary>
class CpuPercentTests : public testing::TestWithParam<CpuCase> { };

TEST_P(CpuPercentTests, NormalizesTotalCapacity)
{
	const auto [Cpu, Elapsed, Processors, Expected] = GetParam();
	Snapshot start, end;
	start.cpu_seconds = 5.0;
	end.cpu_seconds   = 5.0 + Cpu;
	end.timestamp     = start.timestamp + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(Elapsed));
	const auto Actual = CpuPercent(start, end, Processors);
	ASSERT_TRUE(Actual.has_value());
	EXPECT_DOUBLE_EQ(*Actual, Expected);
}

INSTANTIATE_TEST_SUITE_P(Capacity, CpuPercentTests, testing::Values(
	                         CpuCase{0.0, 1.0, 8U, 0.0}, CpuCase{1.0, 1.0, 8U, 12.5},
	                         CpuCase{8.0, 1.0, 8U, 100.0}, CpuCase{4.0, 2.0, 8U, 25.0}, CpuCase{9.0, 1.0, 8U, 100.0}),
                         ([](const testing::TestParamInfo<CpuCase>& info) {
	                         constexpr const char* Names[]{"Idle", "OneBusyCore", "AllBusyCores", "PartialLoad", "CounterRounding"};
	                         return Names[info.index];
	                         }));

/// <summary>Entrées empêchant une interprétation fiable des compteurs CPU.</summary>
class InvalidCpuPercentTests : public testing::TestWithParam<int> { };

TEST_P(InvalidCpuPercentTests, ReportsUnavailable)
{
	Snapshot start, end;
	start.cpu_seconds       = 1.0;
	end.cpu_seconds         = 2.0;
	end.timestamp           = start.timestamp + std::chrono::seconds(1);
	unsigned int processors = 8;
	switch (GetParam()) {
		case 0: start.cpu_seconds.reset();
			break;
		case 1: end.cpu_seconds.reset();
			break;
		case 2: processors = 0;
			break;
		case 3: end.timestamp = start.timestamp;
			break;
		case 4: end.timestamp = start.timestamp - std::chrono::seconds(1);
			break;
		case 5: end.cpu_seconds = 0.0;
			break;
		default: FAIL() << "Scénario inconnu";
	}
	EXPECT_FALSE(CpuPercent(start, end, processors).has_value());
}

INSTANTIATE_TEST_SUITE_P(Unavailable, InvalidCpuPercentTests, testing::Range(0, 6),
                         ([](const testing::TestParamInfo<int>& info) {
	                         constexpr const char* Names[]{"MissingStart", "MissingEnd", "UnknownCapacity", "ZeroDuration", "NegativeDuration",
	                         "DecreasingCounter"};
	                         return Names[info.index];
	                         }));

TEST(ResourceMonitoringTests, ReadsNativeCounters)
{
	// Aucun délai ni seuil de performance : seules la disponibilité et la monotonie sont vérifiées.
	const Snapshot Start = Capture();
	const Snapshot End   = Capture();
	ASSERT_TRUE(Start.cpu_seconds.has_value());
	ASSERT_TRUE(End.cpu_seconds.has_value());
	EXPECT_GE(*End.cpu_seconds, *Start.cpu_seconds);
	EXPECT_GE(End.timestamp, Start.timestamp);
	ASSERT_TRUE(End.resident_bytes.has_value());
	EXPECT_GT(*End.resident_bytes, 0U);
	EXPECT_GT(LogicalProcessorCount(), 0U);
}
} // namespace
} // namespace palm_tracer::monitoring
