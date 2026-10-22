#include <cstdint>
#include <mutex>
#include <omp.h>

#include "ExportAPI.hpp"
#include "Palm.hpp"
#include "TrackAnalysis.hpp"
#include "fit/AstigFit.hpp"
#include "localization/ATrous.hpp"
#include "signal_processing/Alignment.hpp"
#include "signal_processing/Matrix.hpp"
#include "tracking/Reconnection.hpp"
#include "tracking/Tracks.hpp"

namespace col = palm_tracer::localization_columns;


static std::once_flag omp_init_flag;

static void InitOpenMP()
{
	omp_set_dynamic(1);
	// éventuellement d’autres réglages :
	//omp_set_num_threads(static_cast<int>(std::thread::hardware_concurrency())); // Fixe un max
	// omp_set_nested(...)
}

//-------------------------------------------------------------------------------------------------
/// <summary>Aligne un stack 3D d'images à l'aide de coefficients de correction fournis.
///
/// Cette fonction applique une correction géométrique (translation / distorsion) sur chaque plan d'un stack volumique, typiquement pour compenser des décalages
/// entre canaux ou entre plans.
/// Les images d'entrée et de sortie sont stockées sous forme de stack 3D : <c>planes</c> plans successifs de taille <c>width × height</c>.
/// </summary>
/// <param name="input">Pointeur vers le buffer d'entrée contenant le stack d'images.
/// Le buffer doit contenir <c>planes × width × height</c> pixels de type <c>uint16_t</c>, stockés plan par plan (plan 0, puis plan 1, etc.).</param>
/// <param name="output">Pointeur vers le buffer de sortie où seront écrites les images alignées.
/// Le buffer doit être pré-alloué et avoir la même taille que <paramref name="input" />.</param>
/// <param name="height">Hauteur (en pixels) de chaque plan 2D.</param>
/// <param name="width">Largeur (en pixels) de chaque plan 2D.</param>
/// <param name="planes">Nombre de plans dans le stack 3D (dimension Z).</param>
/// <param name="factors"> Tableau de coefficients de correction d'alignement.
/// Il doit contenir exactement 20 valeurs (<c>2 × 10</c>), organisées en 2 lignes de 10 coefficients (une ligne par axe) tels qu'attendus par
/// <c>palm_tracer::Alignment</c>.</param>
/// <param name="upsampling">Facteur de sur-échantillonnage utilisé pour le calcul sub-pixel de l'alignement.
/// Une valeur supérieure à 1 permet un ajustement plus fin mais augmente le coût de calcul.</param>
PT_API void PT_CALL Alignment(const uint16_t* input, uint16_t* output, const uint64_t height,
                              const uint64_t width, const uint64_t planes, const double* factors, const uint64_t upsampling)
{
	std::call_once(omp_init_flag, InitOpenMP);

	std::array<std::array<double, 10>, 2> f{};
	uint64_t idx = 0;
	for (uint64_t i = 0; i < 2; ++i) { for (uint64_t j = 0; j < 10; ++j, ++idx) { f[i][j] = factors[idx]; } }

	palm_tracer::Alignment align;
	const palm_tracer::Vec3I Dim = {.x = width, .y = height, .z = planes};
	align.Align(input, output, Dim, f, upsampling);
}

//-------------------------------------------------------------------------------------------------
/// <summary>Calcule automatiquement un seuil d'intensité pour la détection de molécules sur une image.
///
/// Cette fonction estime un seuil global en analysant la distribution d'intensité de l'image et en utilisant les paramètres d'ajustement <paramref
/// name="fit_params" />.
/// Elle est pensée pour être utilisée avant un appel à <c>Localization</c>.</summary>
/// <param name="input">Pointeur vers le buffer d'entrée contenant le stack d'images.
/// Le buffer doit contenir <c>planes × width × height</c> pixels de type <c>uint16_t</c>, stockés plan par plan (plan 0, puis plan 1, etc.).</param>
/// <param name="height">Hauteur (en pixels) de chaque plan 2D.</param>
/// <param name="width">Largeur (en pixels) de chaque plan 2D.</param>
/// <param name="fit_params">Tableau de paramètres de fit utilisés par l'algorithme de seuil automatique.
/// Le contenu exact dépend de l'implémentation de <c>palm_tracer::Palm</c> (par exemple paramètres de PSF, bruit, etc.).</param>
/// <returns>Seuil d'intensité calculé automatiquement, exprimé dans les mêmes unités que les valeurs de <paramref name="input" /> (converties en
/// <c>double</c>).</returns>
PT_API double PT_CALL AutoThreshold(const uint16_t* input, const uint64_t height, const uint64_t width, const double* fit_params)
{
	std::call_once(omp_init_flag, InitOpenMP);

	const auto PointsListSize = static_cast<uint64_t>(0.2 * static_cast<double>(height * width) * col::Count);
	std::vector<double> pointsList(PointsListSize);
	const palm_tracer::Vec3I Dim = {.x = width, .y = height, .z = 1};

	auto palm = palm_tracer::Palm(input, pointsList.data(), PointsListSize, Dim, 0, 10, 0, fit_params);
	return palm.AutoThreshold();
}

//-------------------------------------------------------------------------------------------------
/// <summary>Effectue la localisation des molécules sur un stack d'images (2D ou 3D).
///
/// Cette fonction détecte les maxima locaux au-dessus d'un seuil, applique éventuellement un watershed, puis réalise un ajustement de PSF (fit) pour chaque
/// spot détecté.
/// Les résultats sont écrits dans <paramref name="output" /> sous forme de tableau plat de <c>double</c> (format défini par
/// <c>palm_tracer::Palm</c>).</summary>
/// <param name="input">Pointeur vers le buffer d'entrée contenant le stack d'images.
/// Le buffer doit contenir <c>planes × width × height</c> pixels de type <c>uint16_t</c>, stockés plan par plan (plan 0, puis plan 1, etc.).</param>
/// <param name="output">Pointeur vers le tableau de sortie contenant les localisations détectées.
/// Le tableau doit être pré-alloué et suffisamment grand pour contenir <paramref name="outputSize"/> valeurs.
/// Le format exact (ordre des colonnes : x, y, z, intensité, etc.) est défini par <c>palm_tracer::Palm</c>.</param>
/// <param name="outputSize">Taille (en nombre de <c>double</c>) du tableau <paramref name="output" />.
/// Si le nombre de localisations dépasse la capacité du tableau, l'excédent sera ignoré.</param>
/// <param name="height">Hauteur (en pixels) de chaque plan 2D.</param>
/// <param name="width">Largeur (en pixels) de chaque plan 2D.</param>
/// <param name="planes">Nombre de plans dans le stack 3D (dimension Z).</param>
/// <param name="thresholdVal">Seuil d'intensité utilisé pour la détection des spots. Les pixels en dessous de ce seuil ne sont pas considérés comme
/// candidats.</param>
/// <param name="watershedRatio">Seuil de fusion des bassins : 0 conserve les maxima distincts (Watershed activé dans Python),
/// 10 autorise leur fusion (Watershed désactivé). Une valeur plus élevée favorise la fusion, et non la séparation.</param>
/// <param name="fit">Mode d'ajustement (fit) à utiliser pour chaque spot. La signification exacte des valeurs dépend de l'implémentation de
/// <c>palm_tracer::Palm</c>.</param>
/// <param name="fit_params">Tableau de paramètres de fit (par exemple paramètres de PSF, initialisation, contraintes).</param>
/// <returns>Nombre de valeurs effectivement écrites dans <paramref name="output" />.</returns>
PT_API uint64_t PT_CALL Localization(const uint16_t* input, double* output, const uint64_t outputSize, const uint64_t height,
                                     const uint64_t width, const uint64_t planes, const double thresholdVal, const double watershedRatio,
                                     const uint64_t fit, const double* fit_params)
{
	std::call_once(omp_init_flag, InitOpenMP);

	const palm_tracer::Vec3I Dim = {.x = width, .y = height, .z = planes};
	const auto Palm              = palm_tracer::Palm(input, output, outputSize, Dim, thresholdVal, watershedRatio, fit, fit_params);
	return Palm.Processing();
}


//-------------------------------------------------------------------------------------------------
/// <summary>Construit des trajectoires 2D à partir d'un tableau contigu de localisations simples, puis exporte le couple <c>(Id, Track)</c> de chaque
/// point.</summary>
/// <param name="input">Tableau d'entrée contigu de doubles contenant 4 valeurs par point : <c>Id</c>, <c>Plane</c>, <c>X</c>, <c>Y</c>.</param>
/// <param name="output">Tableau de sortie contigu préalloué en entiers 64 bits. Il reçoit deux valeurs par point : <c>Id</c>, puis identifiant de
/// trajectoire.</param>
/// <param name="nPoints">Nombre de points présents dans <paramref name="input"/>.</param>
/// <param name="maxDistance">Distance maximale autorisée entre deux points de deux plans consécutifs pour créer une liaison.</param>
/// <remarks>Écrit <c>2 × nPoints</c> valeurs dans <paramref name="output"/> ; cet export ne renvoie pas de valeur.</remarks>
PT_API void PT_CALL Tracking(const double* input, uint64_t* output, const uint64_t nPoints, const double maxDistance)
{
	std::call_once(omp_init_flag, InitOpenMP);

	palm_tracer::tracking::Tracks::Tracking(input, output, nPoints, maxDistance);
}


//-------------------------------------------------------------------------------------------------
/// <summary>Réassocie les trajectoires interrompues à cause du blinking (reconnexion des blinkings).
///
/// Cette fonction charge des trajectoires existantes, tente de reconnecter les segments séparés par des périodes de clignotement (blinking), puis écrit le
/// résultat dans <paramref name="output"/>.
/// </summary>
/// <param name="input">Pointeur vers le tableau d'entrée contenant les trajectoires initiales, tel que produit par <c>Tracking</c> ou une étape
/// précédente.</param>
/// <param name="output">Pointeur vers le tableau de sortie contenant les trajectoires après reconnexion des blinkings. Le buffer doit être pré-alloué.</param>
/// <param name="nRow">Nombre de lignes (trajets / points) dans le tableau d'entrée. La signification exacte dépend du format adopté dans
/// <c>Tracks::LoadTracks</c>.</param>
/// <param name="pixelSize">Taille d'un pixel, en unités physiques (par exemple µm/pixel). Utilisée pour convertir les distances en unités physiques lors de la
/// reconnexion.</param>
/// <param name="mode">Mode de diffusion des cibles. La signification des différents modes est définie par <c>palm_tracer::blinking_mode</c>.</param>
/// <param name="maxDuration">Durée maximale (en nombre de frames) autorisée entre deux segments pour être considérés comme appartenant à la même
/// trajectoire.</param>
/// <param name="maxSpeed">Vitesse maximale autorisée des cibles pour reconnecter deux segments séparés par un blinking.</param>
/// <returns>Nombre de doubles écrits dans <paramref name="output"/>, par lignes de <c>tracking_columns::Count</c> valeurs.</returns>
PT_API uint64_t PT_CALL BlinkingReconnection(const double* input, double* output, const uint64_t nRow, const double pixelSize,
                                             const uint64_t mode, const uint64_t maxDuration, const double maxSpeed)
{
	std::call_once(omp_init_flag, InitOpenMP);

	palm_tracer::tracking::Tracks tracks;
	tracks.LoadTracks(input, nRow, pixelSize);
	palm_tracer::tracking::reconnection::Reconnect(tracks.GetTracks(), mode, maxDuration, maxSpeed);
	return tracks.FillOutput(output);
}

//-------------------------------------------------------------------------------------------------
/// <summary>Analyse les trajectoires (MSD, diffusion instantanée, fits).
///
/// Cette fonction charge un ensemble de trajectoires et, en fonction des drapeaux
/// <paramref name="isMSD" />, <paramref name="isInd" /> et <paramref name="fitMode" />,
/// calcule :
/// <list type="bullet">
/// <item><description>les MSD pour chaque trajectoire,</description></item>
/// <item><description>les coefficients de diffusion instantanée,</description></item>
/// <item><description>et/ou des métriques de fit (linéaire, puissance, exponentiel).</description></item>
/// </list>
/// Les résultats sont encodés sous forme de tableaux irréguliers (voir <c>palm_tracer::matrix::SendIrregularArray</c>).</summary>
/// <param name="input">Pointeur vers le tableau d'entrée contenant les trajectoires, au format attendu par <c>Tracks::LoadTracks</c>.</param>
/// <param name="outputMSD">Pointeur vers le tableau de sortie pour les MSD, si <paramref name="isMSD" /> est vrai.</param>
/// <param name="outputInd">Pointeur vers le tableau de sortie pour les coefficients de diffusion instantanée, si <paramref name="isInd" /> est vrai.</param>
/// <param name="outputFit">Pointeur vers le tableau de sortie pour les résultats des fits (selon <paramref name="fitMode" />).</param>
/// <param name="nRow">Nombre de lignes dans le tableau d'entrée (nombre de points ou de trajectoires, selon la convention de <c>LoadTracks</c>).</param>
/// <param name="isMSD">Si vrai, calcule et remplit <paramref name="outputMSD" /> avec les MSD de toutes les trajectoires.</param>
/// <param name="isInd">Si vrai, calcule et remplit <paramref name="outputInd" /> avec les coefficients de diffusion instantanée pour chaque
/// trajectoire.</param>
/// <param name="is3D">Indique si les trajectoires sont 3D (<c>true</c>) ou 2D (<c>false</c>), ce qui affecte le calcul des métriques.</param>
/// <param name="pixelSize">Taille d'un pixel en unités physiques (par exemple µm/pixel).</param>
/// <param name="exposureTime">Temps d'exposition (ou intervalle temporel entre frames), en secondes.</param>
/// <param name="fitMode">Mode de fit à appliquer :
/// <list type="bullet">
/// <item><description><c>0</c> : aucun fit n'est calculé,</description></item>
/// <item><description><c>1</c> : fit linéaire,</description></item>
/// <item><description><c>2</c> : fit de type puissance,</description></item>
/// <item><description><c>3</c> : fit exponentiel.</description></item>
/// </list></param>
/// <param name="fitParams">Tableau de paramètres pour le fit. L'élément <c>fitParams[0]</c> doit contenir la longueur du fit (<c>fitLength</c>, nombre de
/// points utilisés pour le fit).
/// Les autres éléments peuvent être utilisés par des variantes futures.</param>
/// <returns>Toujours <c>true</c> si l'exécution s'est terminée sans exception. En cas d'erreur interne, une exception C++ peut être levée et doit être gérée
/// côté appelant .NET/Python.</returns>
PT_API bool PT_CALL TrackAnalysis(const double* input, double* outputMSD, double* outputInd, double* outputFit, const uint64_t nRow,
                                  const bool isMSD, const bool isInd, const bool is3D,
                                  const double pixelSize, const double exposureTime, const uint64_t fitMode, const double* fitParams)
{
	std::call_once(omp_init_flag, InitOpenMP);
	return palm_tracer::cpu::AnalyzeTracks(input, outputMSD, outputInd, outputFit, nRow,
	                                       isMSD, isInd, is3D, pixelSize, exposureTime, fitMode, fitParams);
}


//-------------------------------------------------------------------------------------------------
/// <summary>Extrait un plan d'ondelettes d'un stack d'images en utilisant une décomposition à trous.
///
/// Pour chaque plan 2D du stack d'entrée, cette fonction calcule une décomposition en ondelettes (algorithme à trous) via
/// <c>palm_tracer::localization::ATrous</c>
/// et copie le niveau demandé dans <paramref name="output" />.</summary>
/// <param name="input">Pointeur vers le buffer d'entrée contenant le stack d'images.
/// Le buffer doit contenir <c>planes × width × height</c> pixels de type <c>uint16_t</c>, stockés plan par plan (plan 0, puis plan 1, etc.).</param>
/// <param name="output">Pointeur vers le buffer de sortie où seront écrites les images alignées.
/// Le buffer doit être pré-alloué et avoir la même taille que <paramref name="input" />.</param>
/// <param name="height">Hauteur (en pixels) de chaque plan 2D.</param>
/// <param name="width">Largeur (en pixels) de chaque plan 2D.</param>
/// <param name="planes">Nombre de plans dans le stack 3D (dimension Z).</param>
/// <param name="levelOutput">Niveau de la décomposition à extraire.
/// La numérotation commence à 1 côté API : <c>levelOutput = 1</c> correspond au premier niveau d'ondelettes.
/// Si la valeur demandée dépasse le nombre de niveaux disponibles, le niveau maximal est utilisé.</param>
PT_API void PT_CALL Wavelett(const uint16_t* input, double* output, const uint64_t height, const uint64_t width, const uint64_t planes, uint64_t levelOutput)
{
	std::call_once(omp_init_flag, InitOpenMP);
	if (planes == 0) { return; }

	if (levelOutput == 0) { levelOutput = 1; } // on force au niveau 1 si l'appelant donne 0
	// if (levelOutput == 0) { throw std::invalid_argument("levelOutput must be >= 1"); } // lever une exception dans une DLL Extern C n'est pas géré
	// correctement par C# ou Python.

	const palm_tracer::Vec2I Dim = {.x = width, .y = height};
	const uint64_t Size          = Dim.GetSize();
	const uint64_t Level         = std::min(levelOutput - 1, static_cast<uint64_t>(palm_tracer::N_PLANES) - 1); // le niveau 1 correspond à l'index 0 du tableau
	const int Threads            = static_cast<int>(std::min(planes, static_cast<uint64_t>(omp_get_max_threads())));

	#pragma omp parallel num_threads(Threads)
	{
		palm_tracer::localization::ATrous decomposition(Dim);
		// Chaque thread conserve son espace de travail entre les plans qui lui sont attribués.
		#pragma omp for schedule(static)
		for (int i = 0; i < static_cast<int>(planes); ++i) {
			const uint64_t PlaneOffset = i * Size;
			const uint16_t* image      = input + PlaneOffset;
			decomposition.Process(image);
			const double* waveletPlane = decomposition.wDecomp_[Level].data();
			std::copy_n(waveletPlane, Size, output + PlaneOffset);
		}
	}
}

//-------------------------------------------------------------------------------------------------
/// <summary>Calcule la calibration d'un modèle d'astigmatisme 3D à partir d'un jeu de points de calibration.</summary>
/// <param name="input">Pointeur vers le tableau d'entrée contenant les points de calibration, au format (\f$\sigma_x, \sigma_y, z\f$).</param>
/// <param name="output">Pointeur vers le buffer de sortie où seront écrites les valeurs du modèle (5 paramètres pour les 2 axes).</param>
/// <param name="size">Nombre de points dans le tableau d'entrée.</param>
/// <param name="pixelSize">Taille d'un pixel en unités physiques (en nanomètre/pixel).</param>
/// <param name="center">Permet de centrer le modèle pour que \f$\sigma_x(0) \approx \sigma_y(0)\f$.</param>
PT_API void PT_CALL Astigmatism3DCalibration(const double* input, double* output, const uint64_t size, const double pixelSize, const bool center = true)
{
	std::call_once(omp_init_flag, InitOpenMP);

	// conversion de l'entrée en vecteur de points 3D.
	const auto N = static_cast<uint64_t>(size);
	std::vector<palm_tracer::Vec3D> dataset(N);
	for (uint64_t i = 0, id = 0; i < N; ++i) {
		dataset[i] = {.x = input[id], .y = input[id + 1], .z = input[id + 2]};
		id += 3;
	}

	// Calculs
	const palm_tracer::fit::AstigFit Fit(pixelSize);
	const palm_tracer::fit::astig_model Model = Fit.Calibrate(dataset, center);

	// Remplissage de la sortie
	for (const auto& Axis : Model) { for (const auto& Param : Axis) { *output++ = Param; } }
}

//-------------------------------------------------------------------------------------------------
/// <summary>Calcule la position axiale Z d'un jeu de points à partir d'un modèle d'astigmatisme 3D.</summary>
/// <param name="input">Paires de sigmas observés (\f$\sigma_x, \sigma_y\f$), en pixels, pour chaque point.</param>
/// <param name="output">Paires [Z en nm, score Q_Z entre 0 et 1] pour chaque point ; 1 indique une correspondance parfaite des sigmas.</param>
/// <param name="size">Nombre de points dans le tableau d'entrée.</param>
/// <param name="pixelSize">Taille d'un pixel en unités physiques (en nanomètre/pixel).</param>
/// <param name="model">Pointeur vers le buffer du modèle (5 paramètres pour les 2 axes).</param>
/// <param name="zMax">Borne absolue Z (nm).</param>
/// <remarks>Le score de compatibilité Q_Z remplace l'ancienne SSE en pixels carrés, sans changer la recherche de Z.
/// Il n'est ni un R2 ni une précision axiale ; voir la documentation de palm_tracer::fit::AstigFit et \ref axial_quality_score.
/// Le format du buffer reste inchangé. La colonne historique Python <c>MSE Z</c> contient désormais ce score, plus grand = meilleur.</remarks>
PT_API void PT_CALL Astigmatism3DEstimation(const double* input, double* output, const uint64_t size, const double pixelSize, const double* model,
                                            const double zMax)
{
	std::call_once(omp_init_flag, InitOpenMP);

	// conversion de l'entrée en vecteur de points 2D.
	const auto N = static_cast<uint64_t>(size);
	std::vector<palm_tracer::Vec2D> datas(N);
	for (uint64_t i = 0, id = 0; i < N; ++i) {
		datas[i] = {.x = input[id], .y = input[id + 1]};
		id += 2;
	}

	const palm_tracer::fit::astig_model Model = {
		std::array{model[0], model[1], model[2], model[3], model[4]},
		std::array{model[5], model[6], model[7], model[8], model[9]}
	};

	const std::array ZLimit = {-zMax, zMax};

	// Calculs
	const palm_tracer::fit::AstigFit Fit(pixelSize);
	const std::vector<palm_tracer::Vec2D> Res = Fit.Estimate(Model, datas, ZLimit);

	// Remplissage de la sortie
	for (const auto& Z : Res) {
		*output++ = Z.x;
		*output++ = Z.y;
	}
}
