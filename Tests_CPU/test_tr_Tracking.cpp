#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <rapidcsv.h>

namespace palm_tracer::tracking {
const std::filesystem::path PATH          = std::filesystem::path(__FILE__).parent_path();
const std::filesystem::path INPUT_PATH    = PATH / "Input";
const std::filesystem::path OUTPUT_PATH   = PATH / "Output";
const std::filesystem::path REF_FILE_PATH = INPUT_PATH / "localizations-parsed.csv";

//-------------------------------------------------------------------------------------------------
inline std::vector<double> ReadVector1D(const std::filesystem::path& path)
{
	// Ouverture du fichier
	std::ifstream file(path);

	if (!file.is_open()) { throw std::runtime_error("Impossible d'ouvrir le fichier : " + path.string()); }

	std::vector<double> values;
	std::string line;

	// Lecture ligne par ligne
	while (std::getline(file, line)) {
		if (line.empty()) { continue; } // Ignore les lignes vides

		values.push_back(std::stod(line));
	}

	return values;
}

//-------------------------------------------------------------------------------------------------
//TEST(TrackingTests, Discontinuous)
//{
//	const std::vector<double> Input = ReadVector1D(REF_FILE_PATH);
//	std::cout << "Input size: " << Input.size() << std::endl;
//
//	std::vector output(688, 0.0); // Pré-allocation du buffer de sortie
//
//	const auto Param = new Parameters();
//
//	Param->maxSpeed_       = 5;
//	Param->costBirthDeath_ = 0.5;
//
//	Tracking tracking(Input.data(), output.data(), 2, 10.0 / 100.0, 64, Param);
//	const size_t Count = tracking.Processing();
//	delete Param;
//}

//-------------------------------------------------------------------------------------------------
} // namespace palm_tracer::tracking
