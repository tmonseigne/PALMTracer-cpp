///-------------------------------------------------------------------------------------------------
/// \file Analysis.cpp
/// <summary>Analyse des trajectoires : MSD, diffusion instantanée et métriques ajustées.</summary>
/// \author Thibaut Monseigne (CNRS).
/// \version 1.3.0.
/// \date 21/10/2026.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///-------------------------------------------------------------------------------------------------
#include "Analysis.hpp"

#include <cmath>
#include <utility>

#include "Fit.hpp"

namespace palm_tracer::tracking::analysis {
//-------------------------------------------------------------------------------------------------
std::vector<double> ComputeMSD(const Track& track, const size_t start, size_t end)
{
	const auto& P = track.points;
	if (end == 0 || end > P.size()) { end = P.size(); }    // Borne de End (pour le comportement par défaut)
	if (end < 2 + start) { return std::vector<double>(); } // Si l'on a moins de 2 éléments à analyser
	const size_t L = P[end - 1].plane_ - P[start].plane_;  // Lag maximum possible
	// if (L ≤ 0) { return std::vector<double>(); } // Impossible : les trajectoires n'ont pas de points sur le même plan et elles sont triées.

	std::vector sum(L, 0.0), res(L, -1.0);
	std::vector<size_t> cnt(L, 0);

	// Accumulation (toutes paires, regroupées par Lag)
	for (size_t i = start; i < end - 1; ++i) {
		const Point& Pi = P[i];
		for (size_t j = i + 1; j < end; ++j) {
			const Point& Pj  = P[j];
			const size_t Lag = (Pj.plane_ - Pi.plane_) - 1;
			if (Lag >= L) { continue; } // Impossible normalement, mais au cas où.
			const double D = Pj.GetDistance(Pi);
			sum[Lag] += D * D;
			cnt[Lag]++;
		}
	}

	// Sortie de taille exacte, Lag vides = -1
	for (size_t i = 0; i < L; ++i) { if (cnt[i] != 0) { res[i] = sum[i] / static_cast<double>(cnt[i]); } }
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> ComputeAllMSD(const std::vector<Track>& tracks)
{
	const auto N = tracks.size(); // Nombre de trajectoires
	std::vector<std::vector<double>> res(N);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const auto Idx     = static_cast<size_t>(i);
		const Track& Track = tracks[Idx];

		if (Track.points.size() < 2) { res[Idx] = std::vector<double>(); }
		else {
			std::vector<double> msd = ComputeMSD(Track);
			msd.insert(msd.begin(), static_cast<double>(Track.id)); // Ajout de l'identifiant de trajectoire en première colonne
			res[Idx] = std::move(msd);
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
size_t GetFitResultSize(const fit_mode mode)
{
	switch (mode) {
		case fit_mode::linear: return 8;       // [L, I, Da, Db, Derr, a, b, err]
		case fit_mode::power: return 9;        // [L, I, Da, Db, Derr, α, a, err, μ(dist)]
		case fit_mode::exponential: return 10; // [L, I, Da, Db, Derr, a, τ, b, err, sqrt(3/4*(b+a))]
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ComputeInstantDiffusion(const Track& track, const ComputesParams& p)
{
	const auto& P = track.points;
	if (p.fitLength < 1 || p.fitLength > P.size()) { return std::vector<double>{}; }
	const size_t NPoint     = p.fitLength + 1; // Il faut prendre un point de plus pour avoir un vecteur de fitLength à la fin du calcul du MSD
	const double SlopeRatio = (p.is3D ? 6.0 : 4.0) * p.exposureTime; // 2*d*Δt
	const auto N            = P.size() - (p.fitLength); // Nombre de points moins la taille de la fenêtre (-1) si Nt = 3, Fit = 2, Np = 3 ⇒ N = 1
	std::vector res(N, -1.0);
	auto fitP = p;
	fitP.mode = fit_mode::linear;

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) {
		Fit fit(fitP);
		auto msd = ComputeMSD(track, i, i + NPoint);               // Il faut prendre un point de plus pour avoir un vecteur de fitLength à la fin
		if (msd.size() > p.fitLength) { msd.resize(p.fitLength); } // En cas de blinking, on supprime les lags plus longs que fitLength
		if (msd.size() < p.fitLength) { res[i] = -1; }             // En cas de fenêtre plus courte, mais impossible géré en amont
		else {
			const double Err = fit.Compute(msd);

			if ((Err < 0) || fit.HasBlink()) { res[i] = -1; }
			// Une erreur négative indique un problème durant le fitting, un blinking "fort" invalide la fenêtre
			else {
				const double S = fit.GetP()[0] / SlopeRatio; // Pente
				res[i]         = S < MinimumInstantDiffusion ? MinimumInstantDiffusion : S;
			}
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ComputeFitMetrics(const Track& track, const ComputesParams& p)
{
	const auto& P = track.points;
	std::vector res(GetFitResultSize(p.mode), -1.0);
	const auto L = static_cast<double>(P.back().plane_ - P.front().plane_);
	size_t idx   = 0;
	res[idx++]   = L + 1; // Delta de longueur de la track (en nombre de plans)
	res[idx]     = 0.0;
	for (const auto& Point : P) { res[idx] += Point.intensity_; } // Calcul de l'intensité totale
	idx++;

	const auto MSD          = ComputeMSD(track);                     // Calcul du MSD
	const double SlopeRatio = (p.is3D ? 6.0 : 4.0) * p.exposureTime; // 2*d*Δt

	// Fit linéaire des N premiers (avec fitLength p[0])
	auto fitP = p;
	fitP.mode = fit_mode::linear;
	Fit fit(fitP);
	double err               = fit.Compute(MSD);
	std::vector<double> pFit = fit.GetP();
	// Une erreur négative indique un problème durant le fitting, un blinking invalide la fenêtre de base
	if ((err < 0) || fit.HasBlink()) { return {}; }

	double s   = pFit[0] / SlopeRatio; // Pente
	res[idx++] = s < MinimumInstantDiffusion ? MinimumInstantDiffusion : s;
	res[idx++] = pFit[1];
	res[idx++] = err;

	// Fit sur l'ensemble des MSD
	fitP.mode      = p.mode;
	fitP.fitLength = MSD.size();

	fit.SetMode(fitP);
	err  = fit.Compute(MSD);
	pFit = fit.GetP();
	if (err < 0) { return {}; } // Si échec, on arrête.
	// (ici le blinking est autorisé, on fait toute la trajectoire, il y en aura forcément en cas de reconnexion)

	// Remplissage suivant le mode
	switch (p.mode) {
		case fit_mode::linear: {
			s          = pFit[0] / SlopeRatio; // Pente
			res[idx++] = s < MinimumInstantDiffusion ? MinimumInstantDiffusion : s;
			res[idx++] = pFit[1];
			res[idx++] = err;

			break;
		}
		case fit_mode::power: {
			res[idx++] = pFit[0]; // Coefficient de puissance
			res[idx++] = pFit[1]; // Coefficient multiplicateur
			res[idx++] = err;

			const auto &Tb = P.back(), Tf = P.front();
			const auto D   = Tb.GetDistance(Tf) / static_cast<double>(Tb.plane_ - Tf.plane_); // Distance moyenne par frame
			res[idx++]     = D / p.exposureTime;                                              // Vitesse moyenne
			break;
		}
		case fit_mode::exponential: {
			if (pFit[1] + pFit[3] <= 0) { return {}; } // Contrainte de positivité.
			res[idx++] = pFit[1];                      // Plateau
			res[idx++] = pFit[2] / p.exposureTime;     // Tau // Divisé par temps d'exposition
			res[idx++] = pFit[3];                      // b
			res[idx++] = err;
			res[idx++] = std::sqrt(0.75 * (pFit[1] + pFit[3]));
			break;
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> ComputeAllInstantDiffusion(const std::vector<Track>& tracks, const ComputesParams& p)
{
	const auto N = tracks.size(); // Nombre de trajectoires
	std::vector<std::vector<double>> res(N);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const auto Idx     = static_cast<size_t>(i);
		const Track& Track = tracks[Idx];

		auto diff = ComputeInstantDiffusion(Track, p);
		diff.insert(diff.begin(), static_cast<double>(Track.id)); // Ajout de l'identifiant de trajectoire en première colonne
		if (diff.size() <= 1) { res[Idx] = std::vector<double>(); }
		else { res[Idx] = std::move(diff); }
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> ComputeAllFitMetrics(const std::vector<Track>& tracks, const ComputesParams& p)
{
	const auto N = tracks.size(); // Nombre de trajectoires
	std::vector<std::vector<double>> res(N);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const auto Idx     = static_cast<size_t>(i);
		const Track& Track = tracks[Idx];

		auto metrics = ComputeFitMetrics(Track, p);
		metrics.insert(metrics.begin(), static_cast<double>(Track.id)); // Ajout de l'identifiant de trajectoire en première colonne
		if (metrics.size() <= 1) { res[Idx] = std::vector<double>(); }
		else { res[Idx] = std::move(metrics); }
	}
	return res;
}
} // namespace palm_tracer::tracking::analysis
