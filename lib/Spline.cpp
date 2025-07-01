#include "Spline.hpp"

#include <array>
#include "Constants.hpp"
#include "palm/Palm.hpp"


namespace palm_tracer {
void Spline::Calibrate(const uint16_t* image, const double* beads, size_t nBeads, size_t height, size_t width, size_t planes, size_t roi, size_t distance, size_t dz, size_t smoothz)
{
	/**
	 Pour faire le fichier de calibration : Les billes doivent être pile sur la hauteur de la Frame (Nframe / 2)
	 les N premiere frame sont au dessus des billes, la frame N/2 est pile sur les billes et les N dernieres frmaes sont en dessous.
	 On bouge d'un pas de N nanomètres à chaque frame. (donc si 51 frame avec un pas de 10 on vas de -250nm à + 250nm)

	========== calibrate3D_g.m ==========
	Calcul des positions des billes (x,y)
	suprresion des billes non valide
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
	// Récupération des positions des billes
}

} // namespace palm_tracer
