#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef INSTANTIATE_TEST_SUITE_P
#define INSTANTIATE_TEST_SUITE_P INSTANTIATE_TEST_CASE_P
#endif

namespace palm_tracer {
using Table = std::vector<std::vector<double>>;

/// <summary>Dossier des donnees partagees avec les tests Python, independant du repertoire courant.</summary>
inline const std::filesystem::path InputPath = std::filesystem::path(__FILE__).parent_path() / "Input";

/// <summary>Lit un CSV numerique avec en-tete et refuse les fichiers absents ou les lignes incompletes.</summary>
inline Table LoadCsv(const std::filesystem::path& path)
{
	std::ifstream stream(path);
	if (!stream) { throw std::runtime_error("Fichier de test introuvable : " + path.string()); }
	std::string line;
	if (!std::getline(stream, line)) { throw std::runtime_error("CSV vide : " + path.string()); }
	const auto Columns = static_cast<size_t>(std::ranges::count(line, ',')) + 1;
	Table result;
	while (std::getline(stream, line)) {
		if (line.empty() || line == "\r") { continue; }
		std::istringstream row(line);
		std::string cell;
		std::vector<double> values;
		while (std::getline(row, cell, ',')) {
			size_t end         = 0;
			const double Value = std::stod(cell, &end);
			if (cell.find_first_not_of(" \t\r", end) != std::string::npos) { throw std::runtime_error("Cellule CSV invalide : " + path.string()); }
			values.push_back(Value);
		}
		if (values.size() != Columns) { throw std::runtime_error("Nombre de colonnes CSV invalide : " + path.string()); }
		result.push_back(std::move(values));
	}
	if (result.empty()) { throw std::runtime_error("CSV sans donnees : " + path.string()); }
	return result;
}

/// <summary>Aplatit les colonnes selectionnees dans l'ordre attendu par l'API C++.</summary>
inline std::vector<double> Flatten(const Table& table, const std::vector<size_t>& columns)
{
	std::vector<double> result;
	result.reserve(table.size() * columns.size());
	for (const auto& Row : table) { for (const size_t Column : columns) { result.push_back(Row.at(Column)); } }
	return result;
}

/// <summary>Compare toutes les dimensions et valeurs d'un tableau, avec diagnostic de ligne et colonne.</summary>
inline void ExpectTableNear(const Table& actual, const Table& expected, const double absolute = 1e-5, const double relative = 0.0)
{
	ASSERT_EQ(actual.size(), expected.size());
	for (size_t row = 0; row < expected.size(); ++row) {
		SCOPED_TRACE("Ligne " + std::to_string(row));
		ASSERT_EQ(actual[row].size(), expected[row].size());
		for (size_t col = 0; col < expected[row].size(); ++col) {
			SCOPED_TRACE("Colonne " + std::to_string(col));
			EXPECT_NEAR(actual[row][col], expected[row][col], absolute + relative * std::abs(expected[row][col]));
		}
	}
}
} // namespace palm_tracer
