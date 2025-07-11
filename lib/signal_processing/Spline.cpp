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
	/**
	 Pour faire le fichier de calibration : Les billes doivent être pile sur la hauteur de la Frame (Nframe / 2)
	 les N premiere frame sont au dessus des billes, la frame N/2 est pile sur les billes et les N dernieres frmaes sont en dessous.
	 On bouge d'un pas de N nanomètres à chaque frame. (donc si 51 frame avec un pas de 10 on vas de -250nm à + 250nm)

	========== calibrate3D_g.m ==========
	Calcul des positions des billes (x,y)
	supression des billes non valide
	midpoint = moitié du nombre de frame

	========== calibrate3D_g.m ==========
	Récupération des ROI sur les N frames pour les N billes (RoiX, RoiY, Z, N Billes)
	sstack = size(beads(1).stack.image);
	beadsStack = zeros(sstack(1), sstack(2), sstack(3), length(beads)) + NaN;
	allstackst = zeros(sstack(1), sstack(2), sstack(3), length(beads)) + NaN;
	Récupération de la frame du milieu (Z / 2)
	halfstoreframes = round((size(beads(1).stack.image, 3) - 1) / 2);
	Calcul d'une psf moyenne sur X,Y,Z
	mstack = mean(beadsStack, 4, 'omitnan');				% Moyenne sur le nombre de billes (une somme serait suffisante vu la normalisation
	mstack = mstack - min(mstack(:), [], 'omitnan');	% Décalage vers 0
	mstack = mstack / sum(mstack(:), 'omitnan');		% Normalisation

	 *
	 */
	SetMembers(size, roi, distance, dz, smoothZ);

	// Récupération des volumes des billes sur l'image Décallé vers zéro (soustraction du minimum de l'image)
	auto beadsStack = ExtractBeads(stack, beads);
	std::cout << beadsStack.size() << "\t" << beadsStack[0].size() << "\n";
	for (size_t i = 0; i < 5; ++i) {
		const size_t RoiIdx = i * beadSize_.x;
		std::cout << "[";
		for (size_t j = 0; j < 5; ++j) { std::cout << beadsStack[0][RoiIdx + j] << "    "; }
		std::cout << "]\n";
	}
	// Alignement et normalisation des PSFs 3D
	std::vector<std::vector<double>> aligned;
}

//-------------------------------------------------------------------------------------------------
void Spline::SetMembers(const Coord3D& size, const size_t roi, const size_t distance, const size_t dz, const size_t smoothZ)
{
	const size_t NewRoi = roi % 2 == 0 ? roi - 1 : roi; // ROI doit-être impair. Pour éviter les effets de bords, on soustrait en cas de nombre pair.

	size_      = size;
	imSize_    = size_.GetPlaneSize();
	stackSize_ = size_.GetSize();

	beadSize_      = Coord3D{.x = NewRoi, .y = NewRoi, .z = size_.z % 2 == 0 ? size_.z - 1 : size_.z}; // Garanti un z impair
	beadImSize_    = beadSize_.GetPlaneSize();
	beadStackSize_ = beadSize_.GetSize();
	midBead_       = beadSize_.GetMid();

	distanceHalf_ = distance / 2;
	dz_           = dz;
	smoothZ_      = smoothZ;

	// Calcul de la zone d'interet centrale de largeur maximum 13 pour la PSF (au-delà, les éléments sont éclipsés de certains calculs)
	constexpr size_t MaxWindowSize = 13;
	const size_t HalfWindow        = std::min(MaxWindowSize / 2, midBead_.x), // division entière
	             StartX            = midBead_.x - HalfWindow,                 // Avec la division entière Start > 0
	             EndX              = midBead_.x + HalfWindow;                 // Avec la division entière End < roi_

	// Génère la plage d’indices du centre de la ROI de taille maximum 13.
	smallXIndices_.clear();              // Remise à 0 du vecteur (en cas d'appels répétés)
	smallXIndices_.reserve(beadSize_.x); // taille max du range intégralité de la psf récupéré
	for (size_t i = StartX; i <= EndX; ++i) { smallXIndices_.push_back(i); }

	// Calcul de la zone d'interet centrale sur Z.
	const size_t StartZ = std::max(0, static_cast<int>(midBead_.z - distanceHalf_)),
	             EndZ   = std::min(static_cast<int>(midBead_.z + distanceHalf_), static_cast<int>(beadSize_.z) - 1);
	// Génère la plage d’indices des plans autour du centre.
	smallZIndices_.clear();              // Remise à 0 du vecteur (en cas d'appels répétés)
	smallZIndices_.reserve(beadSize_.z); // taille max du range intégralité de la psf récupéré
	for (size_t i = StartZ; i <= EndZ; ++i) { smallZIndices_.push_back(i); }

	//Dimansions de la zone d'interet centrale.
	smallBeadSize_      = Coord3D{.x = smallXIndices_.size(), .y = smallXIndices_.size(), .z = smallZIndices_.size()};
	smallBeadImSize_    = smallBeadSize_.GetPlaneSize();
	smallBeadStackSize_ = smallBeadSize_.GetSize();
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Spline::ExtractBeads(const uint16_t* stack, const std::vector<Coord>& beads) const
{
	const size_t Nbeads   = beads.size();
	const uint16_t Minima = *std::min_element(stack, stack + (stackSize_));
	std::vector beadsStack(Nbeads, std::vector(beadStackSize_, 0.0));

	// Extraction des PSF individuelles
	#pragma omp parallel for
	for (int b = 0; b < Nbeads; ++b) {
		const auto& [Xb, Yb] = beads[b];
		auto& s              = beadsStack[b];

		Coord3D posBead{};
		for (posBead.z = 0; posBead.z < beadSize_.z; ++posBead.z) {
			Coord3D pos{};
			pos.z = posBead.z;
			for (posBead.y = 0; posBead.y < beadSize_.y; ++posBead.y) {
				pos.y = Yb - midBead_.y + posBead.y;
				for (posBead.x = 0; posBead.x < beadSize_.x; ++posBead.x) {
					pos.x                          = Xb - midBead_.x + posBead.x;
					s[posBead.GetIndex(beadSize_)] = static_cast<double>(stack[pos.GetIndex(size_)] - Minima);
				}
			}
		}
	}

	// Calcul de la PSF moyenne (somme avec NaN ignorés)
	std::vector mu(beadStackSize_, 0.0);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(beadStackSize_); ++i) {
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
		for (size_t i = 0; i < beadStackSize_; ++i) {
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
std::vector<double> Spline::GetSubVolume(const std::vector<double>& volume, const bool smallZ) const
{
	std::vector<double> res;
	if (smallZ) {
		res.reserve(smallBeadStackSize_); // On réserve l'espace.
		for (const size_t Z : smallZIndices_) {
			const size_t ZIdx = Z * beadImSize_;
			for (const size_t Y : smallXIndices_) {
				const size_t YIdx = ZIdx + Y * beadSize_.x;
				for (const size_t X : smallXIndices_) { res.push_back(volume[YIdx + X]); }
			}
		}
	}
	else {
		res.reserve(smallBeadImSize_ * beadSize_.z); // On réserve l'espace.
		for (size_t z = 0; z < beadSize_.z; ++z) {
			const size_t ZIdx = z * beadImSize_;
			for (const size_t Y : smallXIndices_) {
				const size_t YIdx = ZIdx + Y * beadSize_.x;
				for (const size_t X : smallXIndices_) { res.push_back(volume[YIdx + X]); }
			}
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
void Spline::RegisterPSF(std::vector<std::vector<double>>& beads)
{
	// On prend les N meilleures billes pour faire une moyenne, mais pas moins de 5 sauf si moins de 5 billes au départ.
	const size_t Nref = std::min<size_t>(std::max<size_t>(5, beads.size() / 2), beads.size());

	//Création de la PSF de référence moyenne sur le volume complet
	std::vector mu(beadStackSize_, 0.0);
	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(beadStackSize_); ++i) {
		double sum = 0.0;
		for (size_t b = 0; b < Nref; ++b) { sum += beads[b][i]; }
		mu[i] = sum / static_cast<double>(Nref);
	}

	// Première passe sur l'intégralité du Z.
	// Alignement
	AlignToRef(mu, beads, false);
	//Normalization chelou

	// Seconde passe sur le sous ensemble d'intérêt au niveau du Z.


	std::vector<Point3D> shifts(beads.size());
	for (size_t i = 0; i < beads.size(); ++i) {
		//get3DShift(mu, beads[i], shifts[i]);
		// mais si on a plus besoin des autres frames on peut couper les Z en dehors ou quand meme utile durant les calculs ?
		// new_b = interp3(b, Xq - shift(k, 2), Yq - shift(k, 1), Zq - shift(k, 3), 'cubic', 0);
		// Si ce n'est pas génant voir pour modifier en place b, peut etre que la transformation des ID sera plus relou
		// que de faire une copie puis ecraser b apres genre b = new_b;
	}
}

//-------------------------------------------------------------------------------------------------
std::vector<Vec3> Spline::AlignToRef(const std::vector<double>& ref, std::vector<std::vector<double>>& beads, const bool smallZ)
{
	std::vector<Vec3> shifts(beads.size());
	const auto SubRef       = GetSubVolume(ref, smallZ); // Stocker la sous-image de la Référence
	const Coord3D SmallSize = {.x = smallBeadSize_.x, .y = smallBeadSize_.y, .z = smallZ ? smallBeadSize_.z : beadSize_.z};

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(beads.size()); ++i) { // Pour chaque bille
		Interpolation3D interp(interpolation_mode::trilinear, true);
		const auto SubB = GetSubVolume(beads[i], smallZ);               // Stocker la sous-image
		shifts[i]       = interp.GetShift(SubRef, SubB, SmallSize);     // Lancer getShift entre la sous ref et la sous image
		beads[i]        = interp.Shift(beads[i], beadSize_, shifts[i]); // Lancer Shift sur l'image de base et Remplacer la bille par ce shift
	}
	return shifts;
}

//-------------------------------------------------------------------------------------------------
void Spline::GetOverlap(const std::vector<std::vector<double>>& beads, const std::vector<Vec3>& shift, std::vector<bool>& good, std::vector<double>& weight,
                        double& norm, const bool smallZ) const
{
	const size_t Nbeads         = beads.size();
	const Coord3D SmallSize     = {.x = smallBeadSize_.x, .y = smallBeadSize_.y, .z = smallZ ? smallBeadSize_.z : beadSize_.z};
	const size_t SmallStackSize = SmallSize.GetSize();

	// Pour un sous volume : const auto Sub = getSubVolume(volume, smallZ);
	// Récupération des sous volumes des billes
	std::vector<std::vector<double>> subBeads(Nbeads);
	for (size_t i = 0; i < Nbeads; ++i) { subBeads[i] = GetSubVolume(beads[i], smallZ); }

	// Moyenne pondérée des sous-volumes valides pour créer une image de référence (refimn)
	std::vector ref(SmallStackSize, 0.0);
	double totalWeight = 0.0;
	for (size_t i = 0; i < Nbeads; ++i) {
		if (!good[i]) { continue; }
		for (size_t j = 0; j < SmallStackSize; ++j) { ref[j] += weight[i] * subBeads[i][j]; }
		totalWeight += weight[i];
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
		if (good[i] && cc[i] > 0.0 && std::isfinite(rms[i])) {
			rms[i] /= cc[i];
			maskValid[i] = true;
		}
	}
	const double a = matrix::MeanWithMask(rms, maskValid);
	double b       = std::numeric_limits<double>::quiet_NaN();
	// Écart-type sur les éléments valides (NaN ignorés)
	if (matrix::CountValid(rms) > 1) {
		double sum_sq = 0.0;
		size_t count  = 0;
		for (size_t i = 0; i < Nbeads; ++i) {
			if (maskValid[i]) {
				const double Diff = rms[i] - a;
				sum_sq += Diff * Diff;
				++count;
			}
		}
		b = (count > 1) ? std::sqrt(sum_sq / static_cast<double>(count - 1)) : std::numeric_limits<double>::quiet_NaN();
	}

	// Si b est invalide (pas assez de points), fallback complet
	const double CutOff = std::isfinite(b) ? a + 2.0 * b : std::numeric_limits<double>::infinity();

	// Maj des bons indices et poids finaux
	for (size_t i = 0; i < Nbeads; ++i) {
		if (std::isnan(rms[i]) || rms[i] > CutOff) {
			good[i]   = false;
			weight[i] = std::numeric_limits<double>::quiet_NaN();
		}
		else { weight[i] = 1.0 / rms[i]; }
	}
}

/*
 *
 *
// Normalise un vecteur (divise par sa moyenne)
void normalizeVolume(std::vector<double>& vol) {
    double mean = 0;
    size_t count = vol.size();
    for (double v : vol) mean += v;
    mean /= static_cast<double>(count);
    if (mean < 1e-15) mean = 1.0; // éviter division par 0
    for (double& v : vol) v /= mean;
}

// Recherche du max dans volume 3D (simple max)
size_t argmax3D(const std::vector<double>& vol) {
    return std::distance(vol.begin(), std::max_element(vol.begin(), vol.end()));
}

// Convertit un index 1D en coord 3D (z,y,x)
Eigen::Vector3i ind2sub3D(size_t ind, size_t Z, size_t Y, size_t X) {
    int x = (int)(ind % X);
    ind /= X;
    int y = (int)(ind % Y);
    ind /= Y;
    int z = (int)(ind);
    return Eigen::Vector3i(z, y, x);
}

// TODO : implémenter interpolation cubic subpixel locale max refinement (getmaxInterp)
// Pour l'instant, juste retourne la coord du max

Vec3 getmaxInterp(const std::vector<double>& V, const Eigen::Vector3d& pos, double dx, int w, size_t Z, size_t Y, size_t X) {
    // À implémenter : interpolation cubique 3D autour de pos
    // Pour l'instant retourne pos sans raffinements
    return Vec3{pos(0), pos(1), pos(2)};
}

// Fonction principale get3Dcorrshift
Vec3 get3Dcorrshift(const std::vector<double>& refim,
                    const std::vector<double>& targetim,
                    size_t roi, size_t zN) {

    // refim et targetim ont taille roi*roi*zN (z,y,x)
    const size_t volsize = roi * roi * zN;

    // Normaliser les volumes
    std::vector<double> ref = refim;
    std::vector<double> target = targetim;
    normalizeVolume(ref);
    normalizeVolume(target);

    // FFT 3D des deux volumes
    std::vector<std::complex<double>> fft_ref, fft_target;
    fft3D(ref, fft_ref, zN, roi, roi);
    fft3D(target, fft_target, zN, roi, roi);

    // Produit de Fourier croisé
    std::vector<std::complex<double>> G;
    crossPowerSpectrum(fft_ref, fft_target, G);

    // iFFT 3D de G
    std::vector<double> Gs;
    ifft3D(G, Gs, zN, roi, roi);

    // Chercher max de Gs (corrélation croisée spatiale)
    size_t maxind = argmax3D(Gs);
    Eigen::Vector3i maxPos = ind2sub3D(maxind, zN, roi, roi);

    // Raffinement subpixel (pas encore implémenté)
    Vec3 shift = getmaxInterp(Gs, maxPos.cast<double>(), 0.05, 1, zN, roi, roi);

    // Calcul du décalage final relatif au centre
    shift.dx = shift.dx / 4.0 - (double)roi / 2.0 + 0.25; // l'offset 1/4 est présent dans MATLAB
    shift.dy = shift.dy / 4.0 - (double)roi / 2.0 + 0.25;
    shift.dz = shift.dz / 4.0 - (double)zN / 2.0 + 0.25;

    return shift;
}
 */
} // namespace palm_tracer
