// Tracking.cpp : Defines the exported functions for the DLL application.
#include "Tracking.h"


void autoCompleteFilename(wchar_t* filename, wchar_t* filenameFolder)
{
	wchar_t buf[1024];

	if (*filename == L'\0') { return; }
	//test if the folder in embedded in the filename
	long i;
	for (i = 0; filename[i] != L'\0' && filename[i] != L':' && filename[i] != L'\\'; i++) { ; }
	if (filename[i] == L'\0') { // no folder in embedded in the filename

		// test if filenameFolder contain a file
		long lastSlash = -1;
		long lastPoint = -1;
		for (long j = 0; filenameFolder[j] != L'\0'; j++) {
			if (filenameFolder[j] == L'\\') {
				lastSlash = j;
				lastPoint = -1;
			}
			if (filenameFolder[j] == L'.') { lastPoint = j; }
		}
		if (lastSlash == -1) { // nothing to do, filename is up-to-date
			return;
		}
		else {
			if (lastPoint == -1) { swprintf(buf, L"%s\\%s", filenameFolder, filename); }
			else {
				for (long k = 0; k <= lastSlash; k++) { buf[k] = filenameFolder[k]; }
				swprintf(buf + lastSlash + 1, L"%s", filename);
			}
			swprintf(filename, L"%s", buf);
		}
	}
}

_declspec(dllexport) long _stdcall Tracking(double* pointsList, unsigned long pointsListSize, double* trc, unsigned long trcSize, double maxSpeed, double dZdX, double minLife, double decrease, unsigned long allowFusiondisso, double costBirth, unsigned long nbDim, unsigned long model, unsigned long nbTimes)
{
	TabTimes* tabTimes;
	int ret;

	PublicParameters* publicParameters = new PublicParameters();

	publicParameters->NbIterations = 300; //By default value of the number of ietrations

	publicParameters->NbSegParams = 13;

	publicParameters->folderOut[0] = L'\0';
	publicParameters->statFile[0] = L'\0';
	publicParameters->fileResASCII[0] = L'\0';

	if (publicParameters->fileResASCII[0] != L'\0') { autoCompleteFilename(publicParameters->fileResASCII, publicParameters->folderOut); }

	publicParameters->nbDim = 2;
	publicParameters->pointsList = pointsList;

	publicParameters->pointsListSize = pointsListSize;
	publicParameters->trc = trc;
	publicParameters->trcSize = trcSize;
	publicParameters->distMaxPastAndFuture = maxSpeed;
	publicParameters->dZDX = dZdX;
	publicParameters->nbMinSuccessiveStructures = minLife;
	publicParameters->temperatureDecrease = decrease;
	publicParameters->allowFusionAndFission = allowFusiondisso;
	publicParameters->costBirthDeath = costBirth;
	publicParameters->model = model;
	publicParameters->likenessIntensityMin = 0.2;

	publicParameters->distMaxNeighbour = 20;
	publicParameters->nbTemporalWindows = 1;
	publicParameters->drawSpeedVector2D = false;
	publicParameters->segImagesInCount = 0;
	publicParameters->nbTimes = nbTimes;
	publicParameters->nbTraj = 0;

	switch (publicParameters->model) {
	case 1: tabTimes = new TabTimesCloserNeighbors(publicParameters);
		break;
	case 2: tabTimes = new TabTimesDistanceIntensity(publicParameters);
		break;
	default: return T_ERR_BAD_MODEL;
	}

	if ((ret = tabTimes->initialisation()) < 0) { return ret; }
	if ((ret = tabTimes->runModel()) < 0) { return ret; }
	if ((ret = tabTimes->giveResult()) < 0) { return ret; }

	delete tabTimes;

	int nb = publicParameters->nbTraj;
	delete publicParameters;

	return nb;
}

_declspec(dllexport) long _stdcall Process(double* pointsList, unsigned long pointsListSize, double* trc, unsigned long trcSize, double maxSpeed, double dZdX, double minLife, double decrease, unsigned long allowFusiondisso, double costBirth, unsigned long nbDim, unsigned long model, unsigned long nbTimes)
{
	TabTimes* tabTimes;
	int ret;

	PublicParameters* publicParameters = new PublicParameters();

	publicParameters->NbIterations = 300; //By default value of the number of ietrations

	publicParameters->NbSegParams = 13;

	publicParameters->folderOut[0] = L'\0';
	publicParameters->statFile[0] = L'\0';
	publicParameters->fileResASCII[0] = L'\0';

	if (publicParameters->fileResASCII[0] != L'\0') { autoCompleteFilename(publicParameters->fileResASCII, publicParameters->folderOut); }

	publicParameters->nbDim = 2;
	publicParameters->pointsList = pointsList;

	publicParameters->pointsListSize = pointsListSize;
	publicParameters->trc = trc;
	publicParameters->trcSize = trcSize;
	publicParameters->distMaxPastAndFuture = maxSpeed;
	publicParameters->dZDX = dZdX;
	publicParameters->nbMinSuccessiveStructures = minLife;
	publicParameters->temperatureDecrease = decrease;
	publicParameters->allowFusionAndFission = allowFusiondisso;
	publicParameters->costBirthDeath = costBirth;
	publicParameters->model = model;
	publicParameters->likenessIntensityMin = 0.2;

	publicParameters->distMaxNeighbour = 20;
	publicParameters->nbTemporalWindows = 1;
	publicParameters->drawSpeedVector2D = false;
	publicParameters->segImagesInCount = 0;
	publicParameters->nbTimes = nbTimes;
	publicParameters->nbTraj = 0;

	switch (publicParameters->model) {
		case 1: tabTimes = new TabTimesCloserNeighbors(publicParameters);
			break;
		case 2: tabTimes = new TabTimesDistanceIntensity(publicParameters);
			break;
		default: return T_ERR_BAD_MODEL;
	}

	if ((ret = tabTimes->initialisation()) < 0) { return ret; }
	if ((ret = tabTimes->runModel()) < 0) { return ret; }
	if ((ret = tabTimes->giveResult()) < 0) { return ret; }

	delete tabTimes;

	int nb = publicParameters->nbTraj;
	delete publicParameters;

	return nb;
}
