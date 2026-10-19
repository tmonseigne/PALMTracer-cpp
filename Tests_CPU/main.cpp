#include "ResourceMonitoring.hpp"

#include <gtest/gtest.h>

#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>
#include <string>

namespace palm_tracer::monitoring {
/// <summary>Formate une mesure sans dépendre de la locale du poste ; conserve explicitement les lectures indisponibles.</summary>
std::string Format(const std::optional<double> value)
{
	if (!value) { return "unavailable"; }
	std::ostringstream stream;
	stream.imbue(std::locale::classic());
	stream << std::fixed << std::setprecision(6) << *value;
	return stream.str();
}

/// <summary>Convertit la mémoire résidente en Mio sans remplacer une lecture échouée par zéro.</summary>
std::optional<double> MemoryMiB(const Snapshot& snapshot)
{
	if (!snapshot.resident_bytes) { return std::nullopt; }
	return static_cast<double>(*snapshot.resident_bytes) / (1024.0 * 1024.0);
}

/// <summary>Mesure chaque test, fixtures comprises, et publie ses ressources après le second relevé.</summary>
/// <remarks>La mémoire concerne le processus entier et peut inclure des allocations conservées par les tests précédents.
/// Les deux relevés ne mesurent pas le pic mémoire intermédiaire. Aucune écriture n'est faite dans la zone mesurée.</remarks>
class ResourceListener final : public testing::EmptyTestEventListener {
public:
	/// <summary>Ouvre la mesure avant la préparation du test.</summary>
	void OnTestStart(const testing::TestInfo&) override { start_ = Capture(); }

	/// <summary>Ferme la mesure après le nettoyage, puis complète la console et les propriétés XML GoogleTest.</summary>
	void OnTestEnd(const testing::TestInfo&) override
	{
		const Snapshot End = Capture();
		const auto Elapsed = Format(std::chrono::duration<double>(End.timestamp - start_.timestamp).count());
		const auto Cpu     = Format(CpuPercent(start_, End, processors_));
		const auto Before  = Format(MemoryMiB(start_));
		const auto After   = Format(MemoryMiB(End));
		testing::Test::RecordProperty("monitor_elapsed_seconds", Elapsed);
		testing::Test::RecordProperty("monitor_cpu_percent", Cpu);
		testing::Test::RecordProperty("monitor_logical_processors", std::to_string(processors_));
		testing::Test::RecordProperty("monitor_rss_start_mib", Before);
		testing::Test::RecordProperty("monitor_rss_end_mib", After);
		std::cout << "[ MONITOR  ] " << Elapsed << " s | CPU total: " << Cpu << " % | RSS: " << Before << " -> " << After << " MiB\n";
	}

private:
	Snapshot start_;                                    ///< Compteurs au début du test.
	unsigned int processors_ = LogicalProcessorCount(); ///< Capacité totale utilisée pour normaliser la charge CPU.
};
} // namespace palm_tracer::monitoring

/// <summary>Initialise GoogleTest et son suivi natif des ressources, sans modifier les algorithmes testés.</summary>
int main(int argc, char** argv)
{
	testing::InitGoogleTest(&argc, argv);
	testing::UnitTest::GetInstance()->listeners().Append(new palm_tracer::monitoring::ResourceListener);
	return RUN_ALL_TESTS();
}
