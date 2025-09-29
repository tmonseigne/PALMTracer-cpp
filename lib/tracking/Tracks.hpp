///-------------------------------------------------------------------------------------------------
///
/// \file Tracks.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------
#pragma once
#include <vector>

#include "Point.hpp"

namespace palm_tracer {
/// <summary>Mode de dispersion des points pour estimer les distances maximales de recherche après scintillement au fil du temps.</summary>
enum class blinking_mode : uint8_t {
	immobile,  ///< Distance max de recherche Fixe avec le temps.
	diffusive, ///< Distance max de recherche connections possible, augmente en racine carrée du temps.
	linear     ///< Distance max de recherche connections possible, augmente linéairement avec le temps.
};

/// <summary>Classe manipulant les trajectoires.</summary>
class Tracks {
public:
	Tracks()  = default; ///< Constructeur de la classe <c>Tracks</c>.
	~Tracks() = default; ///< Destructeur de la classe <c>Tracks</c>.

	/// <summary>Charge un tableau de trajectoires précalculé.</summary>
	/// <param name="input">Tableau de trajectoires d'entrée au format tableau 2D.</param>
	/// <param name="nRow">Nombre de lignes du tableau (le nombre de colonnes est défini ici : <see cref="N_TRACKING"/>).</param>
	void LoadTracking(const double* input, size_t nRow);

	/// <summary>Remplit un tableau de trajectoiresau format tableau 2D.</summary>
	/// <param name="output">Tableau de trajectoires de sortie.</param>
	size_t FillOutput(double* output) const;

	// TODO à terme le tracking sera ici avec une meilleure méthode que le random
	// void Tracking(const double* points, double* tracks, const size_t minLife, const double decrease, const size_t ntimes, Parameters* params);

	/// <summary>Reconnecte des trajectoires séparées par scintillement.</summary>
	/// <param name="mode">Mode de dispersion des points <see cref="blinking_mode"/>.</param>
	/// <param name="maxDuration">Durée maximale d'un scintillemnt.</param>
	/// <param name="maxSpeed">Vitesse maximale d'un point entre deux plans (en pixel).</param>
	void BlinkingReconnection(size_t mode, size_t maxDuration, double maxSpeed);

	/// <summary> Surcharge de l'operateur <<.</summary>
	/// <param name="os">Flux de sortie.</param>
	/// <param name="tracks">Objet à afficher</param>
	/// <returns>Flux modifié.</returns>
	friend std::ostream& operator<<(std::ostream& os, const Tracks& tracks);

protected:
	using track_t = std::vector<Point>; ///< Alias pour représenter une trajectoire composée de points.
	std::vector<track_t> t_;            ///< Liste des Trajectoires.

	void connect(size_t track1, size_t track2);

	static double getMaxDistance(size_t mode, double duration, double maxSpeed);
};

} // namespace palm_tracer
