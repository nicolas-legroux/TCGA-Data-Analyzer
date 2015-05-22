#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "typedefs.hpp"
#include "correlationMatrix.hpp"
#include "stats.hpp"

using namespace std;

void prepareData(std::vector<std::vector<double>> &data,
		std::vector<SampleIdentifier> &sampleIdentifiers,
		CancerPatientIDList &cancerPatientIDList,
		const PatientList &patientControlList,
		const PatientList &patientTumorList, const RNASeqData &controlData,
		const RNASeqData &tumorData) {

	cout << endl << "Preparing data for correlation computation... " << flush;

	int countPatients = 0;

	for (const auto &kv : tumorData) {
		string cancerName = kv.first;
		cancerPatientIDList.insert(
				make_pair(cancerName + "-" + "Tumor", vector<int>()));
		cancerPatientIDList.insert(
				make_pair(cancerName + "-" + "Control", vector<int>()));

		unsigned int numberOfGenes = kv.second.size();

		for (unsigned int j = 0; j < controlData.at(cancerName).at(0).size();
				++j) {

			cancerPatientIDList[cancerName + "-" + "Control"].push_back(
					countPatients);
			sampleIdentifiers.push_back(
					SampleIdentifier(cancerName, false,
							patientControlList.at(cancerName).at(j)));
			vector<double> patientData(numberOfGenes);
			for (unsigned int k = 0; k < numberOfGenes; ++k) {
				patientData[k] = controlData.at(cancerName).at(k).at(j);
			}
			data.push_back(patientData);
			countPatients++;
		}

		for (unsigned int j = 0; j < tumorData.at(cancerName).at(0).size();
				++j) {

			cancerPatientIDList[cancerName + "-" + "Tumor"].push_back(
					countPatients);
			sampleIdentifiers.push_back(
					SampleIdentifier(cancerName, true,
							patientTumorList.at(cancerName).at(j)));
			vector<double> patientData(numberOfGenes);
			for (unsigned int k = 0; k < numberOfGenes; ++k) {
				patientData[k] = tumorData.at(cancerName).at(k).at(j);
			}
			data.push_back(patientData);
			countPatients++;
		}
	}

	cout << "Done." << endl;
}

std::vector<double> pearson(std::vector<std::vector<double>> &data) {
	cout << endl << "Computing Pearson Correlation... " << flush;
	return computePearsonCorrelation(data);
}

std::vector<double> spearman(std::vector<std::vector<double>> &data) {
	cout << endl << "Computing Spearman Correlation... " << flush;
	return computeSpearmanCorrelation(data);
}

std::vector<double> euclidean(std::vector<std::vector<double>> &data) {
	return computePairwiseEuclideanDistance(data);
}

std::vector<double> manhattan(std::vector<std::vector<double>> &data) {
	return computePairwiseManhattanDistance(data);
}

void exportCorrelationMatrix(const std::vector<double> &correlationMatrix,
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
			matrixOutputStream << correlationMatrix[i * N + j] << "\t";
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

void exportClassStats(const std::vector<double> &correlationMatrix,
		const CancerPatientIDList &cancerPatientIDList,
		const vector<SampleIdentifier> &sampleIdentifiers,
		const string &filemaneCorrelationMeans) {

	cout << endl << "Exporting class stats..." << flush;

	vector<string> classes;
	unsigned int N = (int) sqrt(correlationMatrix.size());

	for (const auto &sampleIdentifer : sampleIdentifiers) {
		string cancerName = sampleIdentifer.cancerName;
		if(sampleIdentifer.isTumor){
			classes.push_back(cancerName + "-Tumor");
		}
		else{
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
					//When I = J : we are comparing the same patients, we know the correlation is 1
					if (I != J) {
						data.push_back(correlationMatrix[N * I + J]);
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
