#include "Spline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>

#include "Interpolation.hpp"
#include "signal_processing/FFT.hpp"
#include "signal_processing/Matrix.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Spline::Calibrate(const uint16_t* stack, const std::vector<Coord>& beads, const Coord3D& size, const size_t roi,
                       const size_t distance, const size_t dz, const size_t smoothZ)
{
	SetMembers(size, roi, distance, dz, smoothZ);

	// Récupération des volumes des billes sur l'image Décallé vers zéro (soustraction du minimum de l'image)
	auto beadsStack = ExtractBeads(stack, beads);
	std::cout << beadsStack.size() << "\t" << beadsStack[0].size() << "\n";
	for (size_t i = 0; i < 5; ++i) {
		const size_t RoiIdx = i * beadDim_.x;
		std::cout << "[";
		for (size_t j = 0; j < 5; ++j) { std::cout << beadsStack[0][RoiIdx + j] << "    "; }
		std::cout << "]\n";
	}
	// Alignement et normalisation des PSFs 3D
	std::vector<std::vector<double>> aligned;
}

//-------------------------------------------------------------------------------------------------
void Spline::SetMembers(const Coord3D& dim, const size_t roi, const size_t distance, const size_t dz, const size_t smoothZ)
{
	const size_t NewRoi = roi % 2 == 0 ? roi - 1 : roi; // ROI doit-être impair. Pour éviter les effets de bords, on soustrait en cas de nombre pair.

	dim_       = dim;
	imSize_    = dim_.GetPlaneSize();
	stackSize_ = dim_.GetSize();

	beadDim_    = Coord3D{.x = NewRoi, .y = NewRoi, .z = dim_.z % 2 == 0 ? dim_.z - 1 : dim_.z}; // Garanti un z impair
	beadImSize_ = beadDim_.GetPlaneSize();
	beadSize_   = beadDim_.GetSize();
	midBead_    = beadDim_.GetMid();

	distanceHalf_ = distance / 2;
	dz_           = dz;
	smoothZ_      = smoothZ;

	// Calcul de la zone d'interet centrale de largeur maximum 13 pour la PSF (au-delà, les éléments sont éclipsés de certains calculs)
	constexpr size_t MaxWindowSize = 13;
	const size_t HalfWindow        = std::min(MaxWindowSize / 2, midBead_.x), // division entière
	             StartX            = midBead_.x - HalfWindow,                 // Avec la division entière Start > 0
	             EndX              = midBead_.x + HalfWindow;                 // Avec la division entière End < roi_

	// Génère la plage d’indices du centre de la ROI de taille maximum 13.
	smallXIndices_.clear();             // Remise à 0 du vecteur (en cas d'appels répétés)
	smallXIndices_.reserve(beadDim_.x); // taille max du range intégralité de la psf récupéré
	for (size_t i = StartX; i <= EndX; ++i) { smallXIndices_.push_back(i); }

	// Calcul de la zone d'interet centrale sur Z.
	const size_t StartZ = std::max(0, static_cast<int>(midBead_.z - distanceHalf_)),
	             EndZ   = std::min(static_cast<int>(midBead_.z + distanceHalf_), static_cast<int>(beadDim_.z) - 1);
	// Génère la plage d’indices des plans autour du centre.
	smallZIndices_.clear();             // Remise à 0 du vecteur (en cas d'appels répétés)
	smallZIndices_.reserve(beadDim_.z); // taille max du range intégralité de la psf récupéré
	for (size_t i = StartZ; i <= EndZ; ++i) { smallZIndices_.push_back(i); }

	zIndices_.resize(beadDim_.z);
	std::iota(zIndices_.begin(), zIndices_.end(), 0); // Remplit avec 0, 1, ..., N-1

	//Dimensions de la zone d'interet centrale.
	smallBeadDim_    = Coord3D{.x = smallXIndices_.size(), .y = smallXIndices_.size(), .z = smallZIndices_.size()};
	smallBeadImSize_ = smallBeadDim_.GetPlaneSize();
	smallBeadSize_   = smallBeadDim_.GetSize();
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Spline::ExtractBeads(const uint16_t* stack, const std::vector<Coord>& beads) const
{
	const size_t Nbeads   = beads.size();
	const uint16_t Minima = *std::min_element(stack, stack + stackSize_);
	std::vector beadsStack(Nbeads, std::vector(beadSize_, 0.0));

	// Extraction des PSF individuelles
	#pragma omp parallel for
	for (int b = 0; b < static_cast<int>(Nbeads); ++b) {
		const auto& [Xb, Yb] = beads[b];
		auto& s              = beadsStack[b];

		Coord3D posBead{};
		for (posBead.z = 0; posBead.z < beadDim_.z; ++posBead.z) {
			Coord3D pos{};
			pos.z = posBead.z;
			for (posBead.y = 0; posBead.y < beadDim_.y; ++posBead.y) {
				pos.y = Yb - midBead_.y + posBead.y;
				for (posBead.x = 0; posBead.x < beadDim_.x; ++posBead.x) {
					pos.x                         = Xb - midBead_.x + posBead.x;
					s[posBead.GetIndex(beadDim_)] = static_cast<double>(stack[pos.GetIndex(dim_)] - Minima);
				}
			}
		}
	}

	// Calcul de la PSF moyenne (somme avec NaN ignorés)
	std::vector mu(beadSize_, 0.0);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(beadSize_); ++i) {
		double sum = 0.0;
		for (size_t b = 0; b < Nbeads; ++b) { sum += beadsStack[b][i]; }
		mu[i] = sum / static_cast<double>(Nbeads);
	}

	// Décalage vers 0 et normalisation de la moyenne
	matrix::Normalisation(mu, true);

	// Calcul des distances quadratiques à la moyenne
	std::vector dstack(Nbeads, 0.0);
	#pragma omp parallel for
	for (int b = 0; b < static_cast<int>(Nbeads); ++b) {
		std::vector<double> tmp = beadsStack[b]; // Copie locale (beadsStack ne doit pas être modifié)

		matrix::Normalisation(tmp, true); // Normalisation individuelle

		// Calcul distance quadratique
		double dist2 = 0.0;
		for (size_t i = 0; i < beadSize_; ++i) {
			const double Dist = tmp[i] - mu[i];
			dist2 += Dist * Dist;
		}
		dstack[b] = dist2;
	}

	// Normalisation des distances (inutile)
	// const double mean_d = std::accumulate(dstack.begin(), dstack.end(), 0.0) / Nbeads;
	// for (auto& d : dstack) { d /= mean_d; }

	// Tri des indices selon la distance à la moyenne
	std::vector<size_t> indices(Nbeads);
	std::iota(indices.begin(), indices.end(), 0);
	std::ranges::sort(indices, [&](const size_t i, const size_t j) { return dstack[i] < dstack[j]; });

	std::vector<std::vector<double>> res(Nbeads);
	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(Nbeads); ++i) { res[i] = std::move(beadsStack[indices[i]]); }
	return res;
}


//-------------------------------------------------------------------------------------------------
std::vector<double> Spline::RegisterPSF(std::vector<std::vector<double>>& beads) const
{
	// On prend les N meilleures billes pour faire une moyenne, mais pas moins de 5 sauf si moins de 5 billes au départ.
	const size_t Nbeads = beads.size();
	const size_t Nref   = std::min<size_t>(std::max<size_t>(5, Nbeads / 2), Nbeads);
	std::vector weights(Nbeads, 1.0); // Poids de chaque bille
	double norm = 0.0;

	std::vector mu(beadSize_, 0.0);
	// Fonciton locale de calcul de moyenne à travers les billes
	auto beadMean = [&](const bool update) {
		#pragma omp parallel for
		for (int i = 0; i < static_cast<int>(beadSize_); ++i) {
			double sum   = 0.0;
			size_t count = 0;

			for (size_t b = 0; b < Nref; ++b) {
				if (std::isnan(weights[b])) { continue; }
				const double Val = beads[b][i];
				if (!std::isnan(Val)) {
					sum += Val;
					++count;
				}
			}
			mu[i] = (count > 0) ? sum / static_cast<double>(count) : (update ? mu[i] : std::numeric_limits<double>::quiet_NaN());
		}
	};

	//Création de la PSF de référence moyenne sur le volume complet
	beadMean(false);

	// Première passe sur l'intégralité du Z.
	const auto Shifts = AlignToRef(mu, beads, false);
	// Alignement (attention, différence du fait de l'interpolation, la suite de l'algorithme est compliqué à comparer)
	const auto BeadsNorm = SoftNormalization(beads);     // Normalization douce
	GetOverlap(BeadsNorm, Shifts, weights, norm, false); // Suppression des outliers

	// Recalcul de la PSF de référence moyenne des billes utilisables sur le volume complet (masi évite d'ajouter des nan
	beadMean(true);

	// Seconde passe sur le sous-ensemble d'intérêt au niveau du Z.
	AlignToRef(mu, beads, true);                                                       // Alignement basé sur le sous-ensemble d'intérêt au niveau du Z.
	beads = SoftNormalization(beads);                                                  // Normalization douce
	std::ranges::fill(weights, 1.0);                                                   // Liste des poids des billes
	for (size_t i = 0; i < 3; ++i) { GetOverlap(beads, Shifts, weights, norm, true); } // Lancement 3 fois pour mettre à jour le poids plusieurs fois
	for (auto& b : beads) { matrix::Scale(b, 1 / norm); }                              // Scale des billes

	beadMean(false);
	return mu;
}

//-------------------------------------------------------------------------------------------------
std::vector<Vec3> Spline::AlignToRef(const std::vector<double>& ref, std::vector<std::vector<double>>& beads, const bool smallZ) const
{
	std::vector<Vec3> shifts(beads.size());
	// Stocker la sous-image de la Référence
	const auto SubRef       = matrix::ExtractSubVolume(ref, beadDim_, smallXIndices_, smallXIndices_, smallZ ? smallZIndices_ : zIndices_);
	const Coord3D SmallSize = {.x = smallBeadDim_.x, .y = smallBeadDim_.y, .z = smallZ ? smallBeadDim_.z : beadDim_.z};

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(beads.size()); ++i) { // Pour chaque bille
		Interpolation3D interp(interpolation_mode::trilinear, true);
		// Stocker la sous-image
		const auto SubB = matrix::ExtractSubVolume(beads[i], beadDim_, smallXIndices_, smallXIndices_, smallZ ? smallZIndices_ : zIndices_);
		shifts[i]       = interp.GetShift(SubRef, SubB, SmallSize);    // Lancer getShift entre la sous ref et la sous image
		beads[i]        = interp.Shift(beads[i], beadDim_, shifts[i]); // Lancer Shift sur l'image de base et Remplacer la bille par ce shift
	}
	return shifts;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Spline::SoftNormalization(const std::vector<std::vector<double>>& beads) const
{
	const size_t CenterSize              = std::min<size_t>(3, midBead_.x);
	std::vector<std::vector<double>> out = beads;

	// Génère la plage d’indices du centre de la ROI de taille 7.
	std::vector<size_t> xr;
	xr.reserve(7);
	for (size_t i = midBead_.x - CenterSize; i <= midBead_.x + CenterSize; ++i) { xr.push_back(i); }

	// 4 iterations de la soft normalisation
	for (size_t iter = 0; iter < 4; ++iter) {
		constexpr double QuantileRatio = 0.75;
		std::vector<double> meanVolume(beadSize_, 0.0);
		std::vector<size_t> validCount(beadSize_, 0);

		// Moyenne des billes (NaN-safe)
		for (const auto& B : out) {
			for (size_t i = 0; i < beadSize_; ++i) {
				if (std::isfinite(B[i])) {
					meanVolume[i] += B[i];
					validCount[i]++;
				}
			}
		}
		for (size_t i = 0; i < beadSize_; ++i) {
			meanVolume[i] = validCount[i] > 0 ? meanVolume[i] / static_cast<double>(validCount[i]) : std::numeric_limits<double>::quiet_NaN();
		}

		auto subMean     = matrix::ExtractSubVolume(meanVolume, beadDim_, xr, xr, smallZIndices_);
		const double Q75 = matrix::Quantile(subMean, QuantileRatio); // Calcul du Quantile

		// Mise à jour des billes
		for (auto& b : out) {
			if (matrix::IsInvalid(b)) { continue; } // S'il ne reste aucune valeur utile
			auto sub = matrix::ExtractSubVolume(b, beadDim_, xr, xr, smallZIndices_);

			// Seuillage par intensité sur le 75e percentile
			std::vector<double> filtered;
			filtered.reserve(sub.size());
			for (size_t i = 0; i < sub.size(); ++i) {
				if (std::isfinite(sub[i]) && !std::isnan(sub[i]) && std::isfinite(subMean[i]) && !std::isnan(subMean[i]) && subMean[i] > Q75) {
					filtered.push_back(sub[i] / subMean[i]);
				}
			}

			const double Factor = matrix::Median(filtered);
			if (Factor > 0) { matrix::Scale(b, 1.0 / Factor); }
		}
	}
	return out;
}

//-------------------------------------------------------------------------------------------------
void Spline::GetOverlap(const std::vector<std::vector<double>>& beads, const std::vector<Vec3>& shift, std::vector<double>& weights,
                        double& norm, const bool smallZ) const
{
	const size_t Nbeads         = beads.size();
	const Coord3D SmallSize     = {.x = smallBeadDim_.x, .y = smallBeadDim_.y, .z = smallZ ? smallBeadDim_.z : beadDim_.z};
	const size_t SmallStackSize = SmallSize.GetSize();

	// Pour un sous volume : const auto Sub = getSubVolume(volume, smallZ);
	// Récupération des sous volumes des billes
	std::vector<std::vector<double>> subBeads(Nbeads);
	for (size_t i = 0; i < Nbeads; ++i) {
		subBeads[i] = matrix::ExtractSubVolume(beads[i], beadDim_, smallXIndices_, smallXIndices_, smallZ ? smallZIndices_ : zIndices_);
	}

	// Moyenne pondérée des sous-volumes valides pour créer une image de référence (refimn)
	std::vector ref(SmallStackSize, 0.0);
	double totalWeight = 0.0;
	for (size_t i = 0; i < Nbeads; ++i) {
		if (std::isnan(weights[i])) { continue; }
		for (size_t j = 0; j < SmallStackSize; ++j) { ref[j] += weights[i] * subBeads[i][j]; }
		totalWeight += weights[i];
	}
	matrix::Scale(ref, 1.0 / totalWeight);

	// Corrélation croisée normalisée
	std::vector cc(Nbeads, 0.0);

	const std::vector<bool> RefMask = matrix::GetValidIndex(ref);
	for (size_t i = 0; i < Nbeads; ++i) {
		std::vector<bool> bMask = matrix::GetValidIndex(subBeads[i]);
		double sumD             = 0.0, sumR = 0.0, sumB = 0.0, valid = 0.0;
		for (size_t j = 0; j < SmallStackSize; ++j) {
			if (!bMask[j] || !RefMask[j]) { continue; }
			sumD += subBeads[i][j] * ref[j];
			sumR += ref[j];
			sumB += subBeads[i][j];
			valid += 1;
		}
		cc[i] = valid * (sumD / (sumR * sumB));
	}

	// Normalisation d’amplitude
	norm = matrix::Max(ref);

	Coord3D cropSize = SmallSize;
	matrix::Crop(ref, cropSize, true, true, false);
	matrix::Scale(ref, 1.0 / norm);

	for (auto& bead : subBeads) {
		cropSize = SmallSize;
		matrix::Crop(bead, cropSize, true, true, false);
		matrix::Scale(bead, 1.0 / norm);
	}

	// RMS de l’erreur (res)
	std::vector rms(Nbeads, 0.0);
	for (size_t i = 0; i < Nbeads; ++i) {
		if (shift[i].x <= 3 && shift[i].y <= 3) { rms[i] = matrix::RmsDifference(ref, subBeads[i]); }
		else { rms[i] = std::numeric_limits<double>::quiet_NaN(); }
	}

	// filtrage robuste
	std::vector maskValid(Nbeads, false);
	for (size_t i = 0; i < Nbeads; ++i) {
		if (!std::isnan(weights[i]) && cc[i] > 0.0 && std::isfinite(rms[i])) {
			rms[i] /= cc[i];
			maskValid[i] = true;
		}
	}
	const double Mu = matrix::MeanWithMask(rms, maskValid);
	double b        = std::numeric_limits<double>::quiet_NaN();
	// Écart-type sur les éléments valides (NaN ignorés)
	if (matrix::CountValid(rms) > 1) {
		double sum_sq = 0.0;
		size_t count  = 0;
		for (size_t i = 0; i < Nbeads; ++i) {
			if (maskValid[i]) {
				const double Diff = rms[i] - Mu;
				sum_sq += Diff * Diff;
				++count;
			}
		}
		b = count > 1 ? std::sqrt(sum_sq / static_cast<double>(count - 1)) : std::numeric_limits<double>::quiet_NaN();
	}

	// Si b est invalide (pas assez de points), fallback complet
	const double CutOff = std::isfinite(b) ? Mu + 2.0 * b : std::numeric_limits<double>::infinity();

	// Maj des bons indices et poids finaux
	for (size_t i = 0; i < Nbeads; ++i) { weights[i] = (std::isnan(rms[i]) || rms[i] > CutOff) ? std::numeric_limits<double>::quiet_NaN() : 1.0 / rms[i]; }
}
} // namespace palm_tracer
