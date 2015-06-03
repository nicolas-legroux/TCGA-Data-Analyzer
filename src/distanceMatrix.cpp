#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cassert>
#include <stdexcept>

#include "distanceMatrix.hpp"
#include "distanceMetrics.hpp"
#include "typedefs.hpp"
#include "stats.hpp"

using namespace std;

string distanceMetricName(const DistanceMetric &distanceMetric) {
	switch (distanceMetric) {
	case DistanceMetric::PEARSON_CORRELATION:
		return "pearson-correlation";
	case DistanceMetric::SPEARMAN_CORRELATION:
		return "spearman-correlation";
	case DistanceMetric::EUCLIDEAN_DISTANCE:
		return "euclidean-distance";
	case DistanceMetric::MANHATTAN_DISTANCE:
		return "manhattan-distance";
	case DistanceMetric::COSINE_SIMILARITY:
		return "cosine-similarity";
	default:
		throw invalid_argument("Unknown distance measure.");
	}
}

std::vector<double> computeDistanceMatrix(
		const std::vector<std::vector<double>> &data, DistanceMetric method) {

	switch (method) {
	case DistanceMetric::PEARSON_CORRELATION:
		return computePairwisePearsonCorrelation(data);
	case DistanceMetric::SPEARMAN_CORRELATION:
		return computePairwiseSpearmanCorrelation(data);
	case DistanceMetric::EUCLIDEAN_DISTANCE:
		return computePairwiseEuclideanDistance(data);
	case DistanceMetric::MANHATTAN_DISTANCE:
		return computePairwiseManhattanDistance(data);
	case DistanceMetric::COSINE_SIMILARITY:
		return computePairwiseCosineSimilarity(data);
	default:
		throw invalid_argument("Unknown distance measure.");
	}
}

void exportDistanceMatrix(const std::vector<double> &distanceMatrix,
		const std::vector<SampleIdentifier> &sampleIdentifiers,
		const std::string &filemaneMatrix,
		const std::string &filenamePatientsIDs,
		const std::string &filenameHeatMapLabels) {

	cout << endl << "Exporting Correlation matrix..." << flush;

	ofstream matrixOutputStream("export/" + filemaneMatrix);
	ofstream patientsOutputStream("export/" + filenamePatientsIDs);

	int N = sampleIdentifiers.size();

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			matrixOutputStream << distanceMatrix[i * N + j] << "\t";
		}
		matrixOutputStream << endl;
	}

	string current = "";
	int countCurrent = 0;
	ofstream outputStreamLabels("export/" + filenameHeatMapLabels);

	for (const SampleIdentifier &sampleIdentifier : sampleIdentifiers) {
		string newCurrent = sampleIdentifier.cancerName + "-"
				+ ((sampleIdentifier.isTumor) ? "Tumor" : "Control");
		if (newCurrent != current) {
			if (countCurrent != 0) {
				outputStreamLabels << current << " " << countCurrent << endl;
			}
			current = newCurrent;
			countCurrent = 1;
		} else {
			countCurrent++;
		}
		patientsOutputStream << sampleIdentifier.toString() << endl;
	}

	outputStreamLabels << current << " " << countCurrent << endl;

	cout << " Done." << endl;
}

void exportClassStats(const std::vector<double> &distanceMatrix,
		const CancerPatientIDList &cancerPatientIDList,
		const vector<SampleIdentifier> &sampleIdentifiers,
		const string &filemaneCorrelationMeans) {

	cout << endl << "Exporting class stats..." << flush;

	vector<string> classes;
	unsigned int N = (int) sqrt(distanceMatrix.size());

	for (const auto &sampleIdentifer : sampleIdentifiers) {
		string cancerName = sampleIdentifer.cancerName;
		if (sampleIdentifer.isTumor) {
			classes.push_back(cancerName + "-Tumor");
		} else {
			classes.push_back(cancerName + "-Control");
		}
	}

	auto end_unique = unique(classes.begin(), classes.end());
	classes.erase(end_unique, classes.end());

	unsigned int n = classes.size();
	vector<double> mean_correlation(n * n);
	vector<double> standard_dev_correlation(n * n);

	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = i; j < n; ++j) {
			vector<double> data;
			for (int I : cancerPatientIDList.at(classes[i])) {
				for (int J : cancerPatientIDList.at(classes[j])) {
					// When I = J : we are comparing the same patients,
                  //  we know the correlation is 1
					if (I != J) {
						data.push_back(distanceMatrix[N * I + J]);
					}
				}
			}

			double mean = computeMean(data);
			double standard_dev = computeStandardDeviation(data);
			mean_correlation[n * i + j] = mean;
			mean_correlation[n * j + i] = mean;
			standard_dev_correlation[n * i + j] = standard_dev;
			standard_dev_correlation[n * j + i] = standard_dev;
		}
	}

	ofstream outputStream("export/" + filemaneCorrelationMeans);
	outputStream << "CLASSES";
	for (const string &s : classes) {
		outputStream << "\t" << s << " (" << cancerPatientIDList.at(s).size()
				<< ")";
	}
	outputStream << endl;

	for (unsigned int i = 0; i < n; ++i) {
		outputStream << classes[i] << " ("
				<< cancerPatientIDList.at(classes[i]).size() << ")";
		for (unsigned int j = 0; j < n; ++j) {
			outputStream << "\t" << mean_correlation[n * i + j] << " ("
					<< standard_dev_correlation[n * i + j] << ")";
		}
		outputStream << endl;
	}

	cout << " Done." << endl;
}
