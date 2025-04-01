/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class TabTimes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/


#include "tab_times_dist_inten.h"

//extern PublicParameters *publicParameters;

/********************************************************************************************************************
Constructuor() of TabTimesDistanceIntensity class.
********************************************************************************************************************/
TabTimesDistanceIntensity::TabTimesDistanceIntensity(PublicParameters* _publicParameters): TabTimes(_publicParameters) { publicParameters = _publicParameters; }

/********************************************************************************************************************
allocTabTimes(int _nbTimes) alloc an array of _nbTimes long 
********************************************************************************************************************/
TimeNodes* TabTimesDistanceIntensity::allocTimeNodes(long _time, long _nbNodes) { return new TimeNodesDistanceIntensity(_time, _nbNodes, publicParameters); }


/********************************************************************************************************************
getLikenessGlobal() compute the totalLikeness of hte system
********************************************************************************************************************/
double TabTimesDistanceIntensity::getLikenessGlobal()
{
	double _likeness = 0.;
	for (int _time = 0; _time < nbTimes; _time++) { _likeness += ((TimeNodesDistanceIntensity*)timeNodes[_time])->getLikenessTime(nbTimes); }

	return _likeness;
}

/********************************************************************************************************************
computeConstantes() compute by convergence the values 
publicParameters->stdIntensity
publicParameters->stdMean
publicParameters->stdDistance
********************************************************************************************************************/
void TabTimesDistanceIntensity::computeConstantes()
{
	publicParameters->stdIntensity = 1000;
	publicParameters->stdMean = 1000;
	publicParameters->stdDistance = 1000;
	double stdMeanComputing = 0;
	double stdDistanceComputing = 0;
	double stdIntensityComputing = 0;


	for (int k = 0; k < 20; k++) {
		long countStdInt = 0;
		long countStdMean = 0;
		long coutStdDist = 0;
		stdMeanComputing = 0;
		stdDistanceComputing = 0;
		stdIntensityComputing = 0;


		for (int _time = 0; _time < nbTimes - 1; _time++) {
			((TimeNodesDistanceIntensity*)timeNodes[_time])->computeConstantSigmas(stdIntensityComputing, countStdInt, stdMeanComputing, countStdMean, stdDistanceComputing, coutStdDist, false); //, true);
		}

		if (stdIntensityComputing > 1e-8 && countStdInt > 0) { publicParameters->stdIntensity = stdIntensityComputing / countStdInt; } //*1.0275;//3.4352;
		publicParameters->stdIntensity = sqrt(publicParameters->stdIntensity);

		if (stdMeanComputing > 1e-8 && countStdMean > 0) { publicParameters->stdMean = stdMeanComputing / countStdMean; } //*1.0275;//3.4352;
		publicParameters->stdMean = sqrt(publicParameters->stdMean);

		if (stdDistanceComputing > 1e-8 && coutStdDist > 0) { publicParameters->stdDistance = stdDistanceComputing / coutStdDist; } //*1.0275;//3.4352;
		publicParameters->stdDistance = sqrt(publicParameters->stdDistance);
	}
	publicParameters->stdDistance /= 3; //attention au /3 c'est just pour augmenter l'influence de la distance
}

/********************************************************************************************************************
updateConstantes() ones the Tracking in compute, compute the constant
********************************************************************************************************************/
void TabTimesDistanceIntensity::updateConstantes()
{
	long cpt = 0;
	double covID = 0;
	double covDM = 0;


	publicParameters->stdDistance = 0;
	publicParameters->stdIntensity = 0;
	publicParameters->stdMean = 0;
	publicParameters->covMeanIntensity = 0;
	for (long t = 0; t < nbTimes; t++) {
		TimeNodes* _localTimeNodes = timeNodes[t];
		for (int n = 0; n < _localTimeNodes->nbNodes; n++) {
			Node* node = _localTimeNodes->tabNodes[n];
			if (node->connectionFuture != NULL) {
				publicParameters->stdDistance += node->distance(node->connectionFuture) * node->distance(node->connectionFuture);
				publicParameters->stdIntensity += (node->info->getIntensity() - node->connectionFuture->info->getIntensity()) * (node->info->getIntensity() - node->connectionFuture->info->getIntensity());
				publicParameters->stdMean += (node->info->getMean() - node->connectionFuture->info->getMean()) * (node->info->getMean() - node->connectionFuture->info->getMean());
				covID += (node->info->getIntensity() - node->connectionFuture->info->getIntensity()) * node->distance(node->connectionFuture);
				publicParameters->covMeanIntensity += (node->info->getIntensity() - node->connectionFuture->info->getIntensity()) * (node->info->getMean() - node->connectionFuture->info->getMean());
				covDM += node->distance(node->connectionFuture) * (node->info->getMean() - node->connectionFuture->info->getMean());
				cpt++;
			}
		}
	}
	publicParameters->stdDistance /= cpt;
	publicParameters->stdIntensity /= cpt;
	publicParameters->stdMean /= cpt;
	publicParameters->covMeanIntensity /= cpt;
	publicParameters->stdIntensity = sqrt(publicParameters->stdIntensity);
	publicParameters->stdDistance = sqrt(publicParameters->stdDistance);
	publicParameters->stdMean = sqrt(publicParameters->stdMean);
}


/********************************************************************************************************************
runModel() is the core of the application, it links the structures according the selected model
********************************************************************************************************************/
int TabTimesDistanceIntensity::runModel()
{
	double error;

	long countLearning = 5;

	for (long iLearning = 0; iLearning < countLearning; iLearning++) {
		if (iLearning == countLearning - 1) { error = 20; }
		else { error = 0; }
		long cpt = 0;

		if (iLearning == 0) {
			computeConstantes();
			for (int _time = 0; _time < nbTimes - 1; _time++) { ((TimeNodesDistanceIntensity*)timeNodes[_time])->makeAllInitializationConnections(); }
		}
		else { updateConstantes(); }

		double totalLikeness = getLikenessGlobal();

		int nbInt = 0;

		int bcl = 0;
		do {
			if (iLearning == countLearning - 1) { error *= 1. - (float)publicParameters->temperatureDecrease / 100.; }

			nbInt++;

			double deltaLikeness = 0;
			for (int i = 0; i < publicParameters->NbIterations; i++) { for (int _time = 0; _time < nbTimes - 1; _time++) { deltaLikeness += ((TimeNodesDistanceIntensity*)timeNodes[_time])->tryRandomChange(error); } }
			totalLikeness += deltaLikeness;

			if (fabs(deltaLikeness) < 1e-8) { cpt++; }
			else { cpt = 0; }

			bcl++;
		}
		while (cpt < 10);
	}

	return RET_OK;
}
