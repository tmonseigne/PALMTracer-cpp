/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class TabTimes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/


#include "TabTimesDistanceIntensity.h"
#include "TimeNodesDistanceIntensity.h"
#include <cstdlib>

//extern Parameters *parameters_;

/********************************************************************************************************************
allocTabTimes(int _nbTimes) alloc an array of _nbTimes long 
********************************************************************************************************************/
TimeNodes* TabTimesDistanceIntensity::allocTimeNodes(const long time, const long nNodes) { return new TimeNodesDistanceIntensity(time, nNodes, parameters_); }


/********************************************************************************************************************
getLikenessGlobal() compute the totalLikeness of hte system
********************************************************************************************************************/
double TabTimesDistanceIntensity::GetLikenessGlobal() const
{
	double res = 0.;
	for (int t = 0; t < nbTimes_; ++t) { res += ((TimeNodesDistanceIntensity*)timeNodes_[t])->GetLikenessTime(nbTimes_); }

	return res;
}

/********************************************************************************************************************
computeConstantes() compute by convergence the values 
parameters_->stdIntensity
parameters_->stdMean
parameters_->stdDistance
********************************************************************************************************************/
void TabTimesDistanceIntensity::computeConstantes() const
{
	parameters_->stdIntensity_ = 1000;
	parameters_->stdMean_ = 1000;
	parameters_->stdDistance_ = 1000;
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


		for (int t = 0; t < nbTimes_ - 1; ++t) {
			((TimeNodesDistanceIntensity*)timeNodes_[t])->ComputeConstantSigmas(stdIntensityComputing, countStdInt, stdMeanComputing, countStdMean,
			                                                                    stdDistanceComputing, coutStdDist, false); //, true);
		}

		if (stdIntensityComputing > 1e-8 && countStdInt > 0) { parameters_->stdIntensity_ = stdIntensityComputing / countStdInt; } //*1.0275;//3.4352;
		parameters_->stdIntensity_ = sqrt(parameters_->stdIntensity_);

		if (stdMeanComputing > 1e-8 && countStdMean > 0) { parameters_->stdMean_ = stdMeanComputing / countStdMean; } //*1.0275;//3.4352;
		parameters_->stdMean_ = sqrt(parameters_->stdMean_);

		if (stdDistanceComputing > 1e-8 && coutStdDist > 0) { parameters_->stdDistance_ = stdDistanceComputing / coutStdDist; } //*1.0275;//3.4352;
		parameters_->stdDistance_ = sqrt(parameters_->stdDistance_);
	}
	parameters_->stdDistance_ /= 3; //attention au /3 c'est just pour augmenter l'influence de la distance
}

/********************************************************************************************************************
updateConstantes() ones the Tracking in compute, compute the constant
********************************************************************************************************************/
void TabTimesDistanceIntensity::updateConstantes() const
{
	long cpt = 0;
	double covID = 0;
	double covDM = 0;


	parameters_->stdDistance_ = 0;
	parameters_->stdIntensity_ = 0;
	parameters_->stdMean_ = 0;
	parameters_->covMeanIntensity_ = 0;
	for (int t = 0; t < nbTimes_; ++t) {
		const TimeNodes* local = timeNodes_[t];
		for (int n = 0; n < local->nNodes_; n++) {
			Node* node = local->nodes_[n];
			if (node->future_ != nullptr) {
				parameters_->stdDistance_ += node->Distance(node->future_) * node->Distance(node->future_);
				parameters_->stdIntensity_ += (node->info_->intensity_ - node->future_->info_->intensity_)
					* (node->info_->intensity_ - node->future_->info_->intensity_);
				parameters_->stdMean_ += (node->info_->mean_ - node->future_->info_->mean_)
					* (node->info_->mean_ - node->future_->info_->mean_);
				covID += (node->info_->intensity_ - node->future_->info_->intensity_) * node->Distance(node->future_);
				parameters_->covMeanIntensity_ += (node->info_->intensity_ - node->future_->info_->intensity_)
					* (node->info_->mean_ - node->future_->info_->mean_);
				covDM += node->Distance(node->future_) * (node->info_->mean_ - node->future_->info_->mean_);
				cpt++;
			}
		}
	}
	parameters_->stdDistance_ /= cpt;
	parameters_->stdIntensity_ /= cpt;
	parameters_->stdMean_ /= cpt;
	parameters_->covMeanIntensity_ /= cpt;
	parameters_->stdIntensity_ = sqrt(parameters_->stdIntensity_);
	parameters_->stdDistance_ = sqrt(parameters_->stdDistance_);
	parameters_->stdMean_ = sqrt(parameters_->stdMean_);
}


/********************************************************************************************************************
runModel() is the core of the application, it links the structures according the selected model
********************************************************************************************************************/
int TabTimesDistanceIntensity::RunModel()
{
	double error;

	constexpr long CountLearning = 5;

	for (long iLearning = 0; iLearning < CountLearning; iLearning++) {
		if (iLearning == CountLearning - 1) { error = 20; }
		else { error = 0; }
		long cpt = 0;

		if (iLearning == 0) {
			computeConstantes();
			for (int t = 0; t < nbTimes_ - 1; ++t) { ((TimeNodesDistanceIntensity*)timeNodes_[t])->MakeAllInitializationConnections(); }
		}
		else { updateConstantes(); }

		double totalLikeness = GetLikenessGlobal();

		int nbInt = 0;

		int bcl = 0;
		do {
			if (iLearning == CountLearning - 1) { error *= 1. - parameters_->decrease_ / 100.; }

			nbInt++;

			double deltaLikeness = 0;
			for (int i = 0; i < Parameters::N_ITERATIONS; i++) {
				for (int t = 0; t < nbTimes_ - 1; ++t) { deltaLikeness += ((TimeNodesDistanceIntensity*)timeNodes_[t])->TryRandomChange(error); }
			}
			totalLikeness += deltaLikeness;

			if (fabs(deltaLikeness) < 1e-8) { cpt++; }
			else { cpt = 0; }

			bcl++;
		}
		while (cpt < 10);
	}

	return RET_OK;
}
