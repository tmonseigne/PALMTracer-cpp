#include "Spline.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>
#include "Constants.hpp"
#include "signal_processing/FFT.hpp"
#include "signal_processing/Matrix.hpp"
#include "palm/Palm.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Spline::Calibrate(const uint16_t* stack, const std::vector<Coord>& beads, const size_t height, const size_t width, const size_t planes,
                       const size_t roi, const size_t distance, const size_t dz, size_t smoothZ)
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
	allstacks = zeros(sstack(1), sstack(2), sstack(3), length(beads)) + NaN;
	allstackst = zeros(sstack(1), sstack(2), sstack(3), length(beads)) + NaN;
	Récupération de la frame du milieu (Z / 2)
	halfstoreframes = round((size(beads(1).stack.image, 3) - 1) / 2);
	Calcul d'une psf moyenne sur X,Y,Z
	mstack = mean(allstacks, 4, 'omitnan');				% Moyenne sur le nombre de billes (une somme serait suffisante vu la normalisation
	mstack = mstack - min(mstack(:), [], 'omitnan');	% Décalage vers 0
	mstack = mstack / sum(mstack(:), 'omitnan');		% Normalisation

	 *
	 */
	setMembers(height, width, planes, roi, distance, dz, smoothZ);

	// Récupération des volumes des billes sur l'image Décallé vers zéro (soustraction du minimum de l'image)
	auto beadsStack = extractBeads(stack, beads);
	std::cout << beadsStack.size() << "\t" << beadsStack[0].size() << "\n";
	for (size_t i = 0; i < 5; ++i) {
		const size_t RoiIdx = i * roi_;
		std::cout << "[";
		for (size_t j = 0; j < 5; ++j) { std::cout << beadsStack[0][RoiIdx + j] << "    "; }
		std::cout << "]\n";
	}
	// Alignement et normalisation des PSFs 3D
	std::vector<std::vector<double>> aligned;
}

//-------------------------------------------------------------------------------------------------
void Spline::setMembers(const size_t height, const size_t width, const size_t planes, const size_t roi,
                        const size_t distance, const size_t dz, const size_t smoothZ)
{
	roi_          = roi % 2 == 0 ? roi - 1 : roi; // roi doit-être impair. Pour éviter les effets de bords, on soustrait en cas de nombre pair.
	roiHalf_      = roi_ / 2;
	roi2_         = roi_ * roi_;
	distance_     = distance;
	distanceHalf_ = distance / 2;
	height_       = height;
	width_        = width;
	planes_       = planes;
	planesHalf_   = planes / 2;
	sizeIm_       = width * height;
	sizeStack_    = sizeIm_ * planes;
	sizeBeads_    = planes * roi2_;
	dz_           = dz;
	smoothZ_      = smoothZ;


	// Calcul de la zone d'interet centrale sur Z.
	zStart_       = std::max(0, static_cast<int>(planesHalf_ - distanceHalf_));
	zEnd_         = std::min(static_cast<int>(planesHalf_ + distanceHalf_), static_cast<int>(planes_) - 1);
	zN_           = zEnd_ - zStart_ + 1;

	// Calcul de la zone d'interet centrale de largeur maximum 13 pour la PSF (au-delà, les éléments sont éclipsés de certains calculs)
	constexpr size_t MaxWindowSize = 13;
	const size_t HalfWindow        = std::min(MaxWindowSize / 2, roiHalf_), // division entière
	             Start             = roiHalf_ - HalfWindow,
	             End               = roiHalf_ + HalfWindow; // Avec la division entière End < roi_

	// Génère la plage d’indices du centre de la ROI de taille maximum 13.
	roiCenterIndices.clear();        // Remise à 0 du vecteur (en cas d'appels répétés)
	roiCenterIndices.reserve(roi2_); // taille max du range intégralité de la psf récupéré
	for (size_t i = Start; i <= End; ++i) {
		const size_t IdxRow = i * roi_;
		for (size_t j = Start; j <= End; ++j) { roiCenterIndices.push_back(IdxRow + j); }
	}
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Spline::extractBeads(const uint16_t* stack, const std::vector<Coord>& beads) const
{
	const size_t Nbeads   = beads.size();
	const uint16_t Minima = *std::min_element(stack, stack + (sizeStack_));
	std::vector allstacks(Nbeads, std::vector(sizeBeads_, 0.0));

	// Extraction des PSF individuelles
	#pragma omp parallel for
	for (int b = 0; b < Nbeads; ++b) {
		const auto& [Xb, Yb] = beads[b];
		auto& stackb         = allstacks[b];

		for (size_t z = 0; z < planes_; ++z) {
			for (size_t dy = 0; dy < roi_; ++dy) {
				for (size_t dx = 0; dx < roi_; ++dx) {
					const size_t X    = Xb - roiHalf_ + dx,
					             Y    = Yb - roiHalf_ + dy,
					             Idx  = (z * height_ + Y) * width_ + X,
					             IdxS = z * roi2_ + dy * roi_ + dx;
					stackb[IdxS]      = static_cast<double>(stack[Idx]) - Minima;
				}
			}
		}
	}

	// Calcul de la PSF moyenne (somme avec NaN ignorés)
	std::vector mu(sizeBeads_, 0.0);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(sizeBeads_); ++i) {
		double sum = 0.0;
		for (size_t b = 0; b < Nbeads; ++b) { sum += allstacks[b][i]; }
		mu[i] = sum / static_cast<double>(Nbeads);
	}

	// Décalage vers 0 et normalisation de la moyenne
	Normalisation(mu, true);

	// Calcul des distances quadratiques à la moyenne
	std::vector dstack(Nbeads, 0.0);
	#pragma omp parallel for
	for (int b = 0; b < static_cast<int>(Nbeads); ++b) {
		std::vector<double> stackb = allstacks[b]; // Copie locale (allstacks ne doit pas être modifié)

		Normalisation(stackb, true); // Normalisation individuelle

		// Calcul distance quadratique
		double dist2 = 0.0;
		for (size_t i = 0; i < sizeBeads_; ++i) {
			const double Dist = stackb[i] - mu[i];
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
	for (int i = 0; i < static_cast<int>(Nbeads); ++i) { res[i] = std::move(allstacks[indices[i]]); }

	return res;
}

//-------------------------------------------------------------------------------------------------
void Spline::registerPSF(std::vector<std::vector<double>>& beads) const
{
	// On prend les N meilleures billes pour faire une moyenne, mais pas moins de 5 sauf si moins de 5 billes au départ.
	const size_t Nref = std::min<size_t>(std::max<size_t>(5, beads.size() / 2), beads.size());

	//Création de la PSF de référence moyenne sur le volume complet
	std::vector mu(sizeBeads_, 0.0);
	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(sizeBeads_); ++i) {
		double sum = 0.0;
		for (size_t b = 0; b < Nref; ++b) { sum += beads[b][i]; }
		mu[i] = sum / static_cast<double>(Nref);
	}

	// Alignement
	// Calcul de la zone d'interet centrale de largeur maximum 13
	constexpr size_t MaxWindowSize = 13;
	const size_t HalfWindow        = std::min(MaxWindowSize, roi_) / 2, // division entière
	             Start             = roiHalf_ - HalfWindow,
	             End               = roiHalf_ + HalfWindow;
	// Génère la plage d’indices [start, end]
	std::vector<size_t> range;
	for (size_t i = Start; i <= End; ++i) { range.push_back(i); }

	std::vector<Point3D> shifts(beads.size());
	for (size_t i = 0; i < beads.size(); ++i) {
		get3DShift(mu, beads[i], shifts[i]);
		// mais si on a plus besoin des autres frames on peut couper les Z en dehors ou quand meme utile durant les calculs ?
		// new_b = interp3(b, Xq - shift(k, 2), Yq - shift(k, 1), Zq - shift(k, 3), 'cubic', 0);
		// Si ce n'est pas génant voir pour modifier en place b, peut etre que la transformation des ID sera plus relou
		// que de faire une copie puis ecraser b apres genre b = new_b;
	}
}

//-------------------------------------------------------------------------------------------------
void Spline::get3DShift(const std::vector<double>& ref, const std::vector<double>& target, Point3D& shift) const
{
	shift = {.x = 0, .y = 0, .z = 0, .v = 0}; // (x, y, z) est le shift et v est la cross-correlation.

	// Création d'image HD (*2)
	constexpr size_t FACTOR = 2;
	const size_t Zn_HD      = FACTOR * zN_, Roi_HD = FACTOR * roi_;
	std::vector<double> ref_HD, target_HD;
	//ref_HD = interp3(ref,FACTOR,'cubic');
	//target_HD = interp3(target,FACTOR,'cubic');
	//ref_HD /= mean(refHD);
	//target_HD /= mean(targetHD);

	// FFT 3D des deux volumes
	std::vector<std::complex<double>> ref_FFT, target_FFT;
	FFT3D(ref_HD, ref_FFT, Zn_HD, Roi_HD, Roi_HD);
	FFT3D(target_HD, target_FFT, Zn_HD, Roi_HD, Roi_HD);

	// Produit de Fourier croisé
	std::vector<std::complex<double>> G;
	CrossPowerSpectrum(ref_FFT, target_FFT, G);

	// iFFT 3D de G
	std::vector<double> Gs;
	IFFT3D(G, Gs, Zn_HD, Roi_HD, Roi_HD);
	const auto MaxG = std::ranges::max_element(Gs);
	shift.v         = *MaxG;
	idToPoint3D(std::distance(Gs.begin(), MaxG), shift, Roi_HD, Roi_HD);

	// Raffinement subpixel (pas encore implémenté)
	// getmaxInterp(Gs, shift, 0.05, 1)

	// Calcul du décalage final relatif au centre (à vérifier, car avec les indices matlab un peu yolo)
	// De plus les calculs ont été fait sur une image HD de Facteur 2 (FACTOR)
	shift.x = (shift.x - static_cast<double>(Roi_HD) + 1) / (2 * FACTOR);
	shift.y = (shift.y - static_cast<double>(Roi_HD) + 1) / (2 * FACTOR);
	shift.z = (shift.z - static_cast<double>(Zn_HD) + 1) / (2 * FACTOR);
}

//-------------------------------------------------------------------------------------------------
void Spline::idToPoint3D(size_t id, Point3D& pos, const size_t x, const size_t y)
{
	pos.x = static_cast<double>(id % x);
	id /= x;
	pos.y = static_cast<double>(id % y);
	id /= y;
	pos.z = static_cast<double>(id);
}

/*
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
