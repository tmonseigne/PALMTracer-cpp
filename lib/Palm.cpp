#include "Palm.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <omp.h>
#include <ranges>
#include <vector>

#include "localization/FrameProcessor.hpp"

namespace palm_tracer {
namespace col = localization_columns;

namespace {
/// <summary>Copie les frames dans l'ordre et attribue les identifiants globaux et locaux.</summary>
/// <param name="frames">Localisations retenues par frame, dans l'ordre de la pile.</param>
/// <param name="output">Tampon capable de contenir toutes les localisations fournies.</param>
/// <returns>Nombre de doubles écrits au format historique.</returns>
int assembleLocalizations(const std::vector<localization::FrameLocalizations>& frames, double* output)
{
	// Copie par plan et attribution de l'identifiant global, du plan et de l'index local.
	size_t count = 0, global = 1;
	for (int i = 0; i < static_cast<int>(frames.size()); ++i) {
		size_t local       = 1;
		const double Plane = i + 1;
		for (const auto& Point : frames[i]) {
			std::ranges::copy(Point, output + count);
			output[count + col::Id]    = static_cast<double>(global);
			output[count + col::Plane] = Plane;
			output[count + col::Index] = static_cast<double>(local);
			++global;
			++local;
			count += col::Count;
		}
	}

	return static_cast<int>(count);
}
} // namespace


//-------------------------------------------------------------------------------------------------
Palm::Palm(const uint16_t* image, double* points, const size_t pointsListSize, const Vec3I& dim,
           const double thresholdVal, const double watershedRatio, const size_t fit, const double* fitParams)
	: stack_(image), points_(points), dim_(dim), size_(dim.GetPlaneSize()),
	  nPointsMaxPerPlane_(dim.z == 0 ? 0 : (pointsListSize / col::Count) / dim.z),
	  thresholdVal_(thresholdVal), watershedRatio_(watershedRatio), fit_(fit), fitParams_(fitParams) { }

//-------------------------------------------------------------------------------------------------
int Palm::Processing() const
{
	if (nPointsMaxPerPlane_ == 0) { return 0; }
	std::vector<localization::FrameLocalizations> allPoints(dim_.z);
	const int Threads = static_cast<int>(std::min(dim_.z, static_cast<size_t>(omp_get_max_threads())));

	#pragma omp parallel num_threads(Threads)
	{
		const Vec2I Dim{.x = dim_.x, .y = dim_.y};
		localization::FrameProcessor processor(Dim, nPointsMaxPerPlane_, watershedRatio_, fit_, fitParams_);
		#pragma omp for
		for (int i = 0; i < static_cast<int>(dim_.z); ++i) { processor.Process(stack_ + (static_cast<size_t>(i) * size_), thresholdVal_, allPoints[i]); }
	}

	return assembleLocalizations(allPoints, points_);
}

//-------------------------------------------------------------------------------------------------
double Palm::AutoThreshold()
{
	localization::FrameLocalizations points;
	std::vector mask(size_, false); //  Création du masque binaire (true = masque actif)
	const double HalfRoi = fitParams_[0] / 2;

	// Calcul de l'écart type initial
	double mu       = std::accumulate(stack_, stack_ + size_, 0.0) / static_cast<double>(size_);
	double variance = 0.0;
	for (size_t i = 0; i < size_; ++i) { variance += (stack_[i] - mu) * (stack_[i] - mu); }
	variance /= static_cast<double>(size_);
	thresholdVal_ = std::sqrt(variance);
	const Vec2I Dim{.x = dim_.x, .y = dim_.y};
	localization::FrameProcessor processor(Dim, nPointsMaxPerPlane_, watershedRatio_, fit_, fitParams_);

	for (size_t i = 0; i < THRESHOLD_ITERATION; ++i) {
		processor.Process(stack_, thresholdVal_, points);
		// Mise à jour du masque
		std::fill(mask.begin(), mask.end(), false);
		for (const auto& P : points) {
			const auto Xmin = static_cast<size_t>(std::max(0, static_cast<int>(P[col::X] - HalfRoi))),
			           Xmax = std::min(dim_.x, static_cast<size_t>(P[col::X] + HalfRoi)),
			           Ymin = static_cast<size_t>(std::max(0, static_cast<int>(P[col::Y] - HalfRoi))),
			           Ymax = std::min(dim_.y, static_cast<size_t>(P[col::Y] + HalfRoi));
			for (size_t y = Ymin; y < Ymax; ++y) { for (size_t x = Xmin; x < Xmax; ++x) { mask[(y * dim_.x) + x] = true; } }
		}
		// Deux passages dans le même ordre, sans recopier les pixels hors masque.
		size_t count = 0;
		double sum   = 0;
		for (size_t j = 0; j < size_; ++j) {
			if (!mask[j]) {
				sum += stack_[j];
				++count;
			}
		}

		// Si nombre de pixels restant insuffisant break
		if (count == 0) { break; }

		// Calcul de l'écart type pour les pixels hors segmentation
		mu       = sum / static_cast<double>(count);
		variance = 0.0;
		for (size_t j = 0; j < size_; ++j) { if (!mask[j]) { variance += (stack_[j] - mu) * (stack_[j] - mu); } }
		variance /= static_cast<double>(count);
		variance = std::sqrt(variance);
		if (std::abs(thresholdVal_ - variance) <= std::numeric_limits<double>::epsilon()) { break; } // Aucun changement
		thresholdVal_ = variance;
	}

	return thresholdVal_;
}

//-------------------------------------------------------------------------------------------------
} // namespace palm_tracer
